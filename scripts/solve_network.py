"""Solve network."""

import pypsa
import re

import numpy as np
import pandas as pd

from pypsa.linopt import get_var, linexpr, define_constraints, join_exprs

from pypsa.linopf import network_lopf, ilopf, join_exprs

from vresutils.benchmark import memory_logger

from helper import override_component_attrs, update_config_with_sector_opts

import logging
logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def add_land_use_constraint(n):

    if 'm' in snakemake.wildcards.clusters:
        _add_land_use_constraint_m(n)
    else:
        _add_land_use_constraint(n)


def _add_land_use_constraint(n):
    #warning: this will miss existing offwind which is not classed AC-DC and has carrier 'offwind'

    for carrier in ['solar', 'onwind', 'offwind-ac', 'offwind-dc']:
        existing = n.generators.loc[n.generators.carrier==carrier,"p_nom"].groupby(n.generators.bus.map(n.buses.location)).sum()
        existing.index += " " + carrier + "-" + snakemake.wildcards.planning_horizons
        n.generators.loc[existing.index,"p_nom_max"] -= existing

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def _add_land_use_constraint_m(n):
    # if generators clustering is lower than network clustering, land_use accounting is at generators clusters

    planning_horizons = snakemake.config["scenario"]["planning_horizons"]
    grouping_years = snakemake.config["existing_capacities"]["grouping_years"]
    current_horizon = snakemake.wildcards.planning_horizons

    for carrier in ['solar', 'onwind', 'offwind-ac', 'offwind-dc']:

        existing = n.generators.loc[n.generators.carrier==carrier,"p_nom"]
        ind = list(set([i.split(sep=" ")[0] + ' ' + i.split(sep=" ")[1] for i in existing.index]))

        previous_years = [
            str(y) for y in
            planning_horizons + grouping_years
            if y < int(snakemake.wildcards.planning_horizons)
        ]

        for p_year in previous_years:
            ind2 = [i for i in ind if i + " " + carrier + "-" + p_year in existing.index]
            sel_current = [i + " " + carrier + "-" + current_horizon for i in ind2]
            sel_p_year = [i + " " + carrier + "-" + p_year for i in ind2]
            n.generators.loc[sel_current, "p_nom_max"] -= existing.loc[sel_p_year].rename(lambda x: x[:-4] + current_horizon)

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def prepare_network(n, solve_opts=None):

    if 'clip_p_max_pu' in solve_opts:
        for df in (n.generators_t.p_max_pu, n.generators_t.p_min_pu, n.storage_units_t.inflow):
            df.where(df>solve_opts['clip_p_max_pu'], other=0., inplace=True)

    if solve_opts.get('load_shedding'):
        n.add("Carrier", "Load")
        n.madd("Generator", n.buses.index, " load",
               bus=n.buses.index,
               carrier='load',
               sign=1e-3, # Adjust sign to measure p and p_nom in kW instead of MW
               marginal_cost=1e2, # Eur/kWh
               # intersect between macroeconomic and surveybased
               # willingness to pay
               # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
               p_nom=1e9 # kW
        )

    if solve_opts.get('noisy_costs'):
        for t in n.iterate_components():
            #if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if 'marginal_cost' in t.df:
                np.random.seed(174)
                t.df['marginal_cost'] += 1e-2 + 2e-3 * (np.random.random(len(t.df)) - 0.5)

        for t in n.iterate_components(['Line', 'Link']):
            np.random.seed(123)
            t.df['capital_cost'] += (1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)) * t.df['length']

    if solve_opts.get('nhours'):
        nhours = solve_opts['nhours']
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760./nhours

    if snakemake.config['foresight'] == 'myopic':
        add_land_use_constraint(n)

    return n


def add_battery_constraints(n):

    chargers_b = n.links.carrier.str.contains("battery charger")
    chargers = n.links.index[chargers_b & n.links.p_nom_extendable]
    dischargers = chargers.str.replace("charger", "discharger")

    if chargers.empty or ('Link', 'p_nom') not in n.variables.index:
        return

    link_p_nom = get_var(n, "Link", "p_nom")

    lhs = linexpr((1,link_p_nom[chargers]),
                  (-n.links.loc[dischargers, "efficiency"].values,
                   link_p_nom[dischargers].values))

    define_constraints(n, lhs, "=", 0, 'Link', 'charger_ratio')


def add_chp_constraints(n):

    electric_bool = (n.links.index.str.contains("urban central")
                     & n.links.index.str.contains("CHP")
                     & n.links.index.str.contains("electric"))
    heat_bool = (n.links.index.str.contains("urban central")
                 & n.links.index.str.contains("CHP")
                 & n.links.index.str.contains("heat"))

    electric = n.links.index[electric_bool]
    heat = n.links.index[heat_bool]

    electric_ext = n.links.index[electric_bool & n.links.p_nom_extendable]
    heat_ext = n.links.index[heat_bool & n.links.p_nom_extendable]

    electric_fix = n.links.index[electric_bool & ~n.links.p_nom_extendable]
    heat_fix = n.links.index[heat_bool & ~n.links.p_nom_extendable]

    link_p = get_var(n, "Link", "p")

    if not electric_ext.empty:

        link_p_nom = get_var(n, "Link", "p_nom")

        #ratio of output heat to electricity set by p_nom_ratio
        lhs = linexpr((n.links.loc[electric_ext, "efficiency"]
                       *n.links.loc[electric_ext, "p_nom_ratio"],
                       link_p_nom[electric_ext]),
                      (-n.links.loc[heat_ext, "efficiency"].values,
                       link_p_nom[heat_ext].values))

        define_constraints(n, lhs, "=", 0, 'chplink', 'fix_p_nom_ratio')

        #top_iso_fuel_line for extendable
        lhs = linexpr((1,link_p[heat_ext]),
                      (1,link_p[electric_ext].values),
                      (-1,link_p_nom[electric_ext].values))

        define_constraints(n, lhs, "<=", 0, 'chplink', 'top_iso_fuel_line_ext')

    if not electric_fix.empty:

        #top_iso_fuel_line for fixed
        lhs = linexpr((1,link_p[heat_fix]),
                      (1,link_p[electric_fix].values))

        rhs = n.links.loc[electric_fix, "p_nom"].values

        define_constraints(n, lhs, "<=", rhs, 'chplink', 'top_iso_fuel_line_fix')

    if not electric.empty:

        #backpressure
        lhs = linexpr((n.links.loc[electric, "c_b"].values
                       *n.links.loc[heat, "efficiency"],
                       link_p[heat]),
                      (-n.links.loc[electric, "efficiency"].values,
                       link_p[electric].values))

        define_constraints(n, lhs, "<=", 0, 'chplink', 'backpressure')

def basename(x):
     return x.split("-2")[0]

def add_pipe_retrofit_constraint(n):
    """Add constraint for retrofitting existing CH4 pipelines to H2 pipelines."""

    gas_pipes_i = n.links.query("carrier == 'gas pipeline' and p_nom_extendable").index
    h2_retrofitted_i = n.links.query("carrier == 'H2 pipeline retrofitted' and p_nom_extendable").index

    if h2_retrofitted_i.empty or gas_pipes_i.empty: return

    link_p_nom = get_var(n, "Link", "p_nom")

    CH4_per_H2 = 1 / n.config["sector"]["H2_retrofit_capacity_per_CH4"]
    fr = "H2 pipeline retrofitted"
    to = "gas pipeline"

    pipe_capacity = n.links.loc[gas_pipes_i, 'p_nom'].rename(basename)

    lhs = linexpr(
        (CH4_per_H2, link_p_nom.loc[h2_retrofitted_i].rename(index=lambda x: x.replace(fr, to))),
        (1, link_p_nom.loc[gas_pipes_i])
    )

    lhs.rename(basename, inplace=True)
    define_constraints(n, lhs, "=", pipe_capacity, 'Link', 'pipe_retrofit')


def add_co2_sequestration_limit(n, sns):

    co2_stores = n.stores.loc[n.stores.carrier=='co2 stored'].index

    if co2_stores.empty or ('Store', 'e') not in n.variables.index:
        return

    vars_final_co2_stored = get_var(n, 'Store', 'e').loc[sns[-1], co2_stores]

    lhs = linexpr((1, vars_final_co2_stored)).sum()

    limit = n.config["sector"].get("co2_sequestration_potential", 200) * 1e6

    opts = snakemake.wildcards.sector_opts.split('-')
    for o in opts:
        if not "seq" in o: continue
        limit = float(o[o.find("seq")+3:]) * 1e6
        break

    name = 'co2_sequestration_limit'
    sense = "<="

    n.add("GlobalConstraint", name, sense=sense, constant=limit,
          type=np.nan, carrier_attribute=np.nan)

    define_constraints(n, lhs, sense, limit, 'GlobalConstraint',
                       'mu', axes=pd.Index([name]), spec=name)

def add_biofuel_constraint(n):

    options = snakemake.wildcards.sector_opts.split('-')
    liquid_biofuel_limit = 0
    biofuel_constraint_type = 'Lt'
    for o in options:
        if "B" in o:
            liquid_biofuel_limit = o[o.find("B") + 1:o.find("B") + 4]
            liquid_biofuel_limit = float(liquid_biofuel_limit.replace("p", "."))
            if len(o) > 7:
                biofuel_constraint_type = o[o.find("B") + 6:o.find("B") + 8]

    print('Liq biofuel minimum constraint: ', liquid_biofuel_limit, ' ', type(liquid_biofuel_limit))

    biofuel_i = n.links.query('carrier == "biomass to liquid"|carrier == "biomass to liquid CC"').index
    biofuel_vars = get_var(n, "Link", "p").loc[:, biofuel_i]
    biofuel_vars_eta = n.links.query('carrier == "biomass to liquid"|carrier == "biomass to liquid CC"').efficiency

    napkership = n.loads.p_set.filter(regex='naphtha for industry|kerosene for aviation|agriculture machinery oil$|shipping oil$').sum() * len(n.snapshots)
    landtrans = n.loads_t.p_set.filter(regex='land transport oil$').sum().sum()

    total_oil_load = napkership+landtrans
    limit = liquid_biofuel_limit * total_oil_load

    lhs = linexpr((biofuel_vars_eta, biofuel_vars)).sum().sum()
    name = 'liquid_biofuel_min'
    sense = '>='
    if biofuel_constraint_type == 'Eq':
        sense = '=='
    elif biofuel_constraint_type == 'Lt':
        sense = '>='

    define_constraints(n, lhs, sense, limit, 'Link', spec=name)


def add_biomass_constraint(n):

    options = snakemake.wildcards.sector_opts.split('-')
    biomass_limit = 0
    for o in options:
        if "biomass" in o:
            biomass_limit = float(o[o.find("biomass") + 7:])# * 1e6 #TWh -> MWh

    hours = list(filter(re.compile(r'^\d+h$', re.IGNORECASE).search, opts))
    hours = [int(s) for s in re.findall(r'\d+',hours[0])]

    print('Adding biomass constraint of max. ', biomass_limit / 1e6, 'TWh')

    biomass_i = n.generators.query('carrier == "solid biomass import" | carrier == "forest residues solid biomass" | carrier == "industry wood residues solid biomass" | carrier == "landscape care solid biomass" | carrier == "manureslurry digestible biomass" | carrier == "sewage sludge digestible biomass" | carrier == "straw digestible biomass"').index
    biomass_vars = get_var(n, "Generator", "p").loc[:, biomass_i]

    limit = biomass_limit / hours[0]

    lhs = linexpr((1,biomass_vars)).sum().sum()
    name = 'biomass_limit'
    sense = '<='

    define_constraints(n, lhs, sense, limit, 'Generator', spec=name)


def add_msw_full_usage(n):
    print('Adding MSW incineration constraint')
    wasteCHP_i = n.links.carrier.filter(regex='waste CHP').index
    wasteCHP_vars = get_var(n, "Link", "p").loc[:, wasteCHP_i]

    lhs = linexpr((1, wasteCHP_vars)).sum().sum()
    rhs = n.generators.p_nom.filter(regex='municipal solid waste').sum() * len(n.snapshots)
    name = 'msw_full_incineration'
    sense = '=='
    define_constraints(n, lhs, sense, rhs, 'Link', spec=name)


def extra_functionality(n, snapshots):
    add_battery_constraints(n)
    add_pipe_retrofit_constraint(n)
    add_co2_sequestration_limit(n, snapshots)

    options = snakemake.wildcards.sector_opts.split('-')
    for o in options:
        if "B" in o:
            print('adding biofuel constraints')
            add_msw_full_usage(n)
            add_biofuel_constraint(n)
        if "biomass" in o:
            print('adding biomass constraint')
            add_biomass_constraint(n)
        if 'CCL' in o:
            print('adding ccl constraints')
            add_ccl_constraints(n)
        if 'convCCL' in o:
            print('adding conventional ccl constraints')
            add_ccl_constraints_conventional(n)


def solve_network(n, config, opts='', **kwargs):
    solver_options = config['solving']['solver'].copy()
    solver_name = solver_options.pop('name')
    cf_solving = config['solving']['options']
    track_iterations = cf_solving.get('track_iterations', False)
    min_iterations = cf_solving.get('min_iterations', 4)
    max_iterations = cf_solving.get('max_iterations', 6)
    keep_shadowprices = cf_solving.get('keep_shadowprices', True)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if cf_solving.get('skip_iterations', False):
        network_lopf(n, solver_name=solver_name, solver_options=solver_options,
                     extra_functionality=extra_functionality,
                     keep_shadowprices=keep_shadowprices,
                     keep_references=True,
                     keep_files=True,
                     **kwargs)
    else:
        ilopf(n, solver_name=solver_name, solver_options=solver_options,
              track_iterations=track_iterations,
              min_iterations=min_iterations,
              max_iterations=max_iterations,
              extra_functionality=extra_functionality,
              keep_shadowprices=keep_shadowprices,
              keep_references=True,
              keep_files=True,
              **kwargs)
    return n


if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'solve_network',
            simpl='',
            opts="",
            clusters="37",
            lv=1.0,
            sector_opts='168H-T-H-B-I-A-solar+p3-dist1',
            planning_horizons="2030",
        )

    logging.basicConfig(filename=snakemake.log.python,
                        level=snakemake.config['logging_level'])

    update_config_with_sector_opts(snakemake.config, snakemake.wildcards.sector_opts)

    tmpdir = snakemake.config['solving'].get('tmpdir')
    if tmpdir is not None:
        from pathlib import Path
        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = snakemake.wildcards.sector_opts.split('-')
    solve_opts = snakemake.config['solving']['options']

    fn = getattr(snakemake.log, 'memory', None)
    with memory_logger(filename=fn, interval=30.) as mem:

        overrides = override_component_attrs(snakemake.input.overrides)
        n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

        n = prepare_network(n, solve_opts)

        n = solve_network(n, config=snakemake.config, opts=opts,
                          solver_dir=tmpdir,
                          solver_logfile=snakemake.log.solver)

        if "lv_limit" in n.global_constraints.index:
            n.line_volume_limit = n.global_constraints.at["lv_limit", "constant"]
            n.line_volume_limit_dual = n.global_constraints.at["lv_limit", "mu"]

        cluster = snakemake.wildcards.clusters
        lv = snakemake.wildcards.lv

        # Save shadow price of biofuel constraint
        #sector_opts = snakemake.wildcards.sector_opts
        #planning_horizon = snakemake.wildcards.planning_horizons[-4:]
        #headerBiofuelConstraint = '{}_lv{}_{}_{}'.format(cluster, lv, sector_opts, planning_horizon)
        #biofuelConstraintFile = snakemake.config['results_dir'] + snakemake.config['run'] + '/biofuelminConstraint.csv'
        #df = pd.DataFrame()
        #print(n.duals["Link"]["df"]["liquid_biofuel_min"])
        #data = pd.DataFrame({f'{headerBiofuelConstraint}': n.duals["Link"]["df"]["liquid_biofuel_min"].values}).T

        #data.to_csv(biofuelConstraintFile, mode='a', header=False)
        ###

        n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
