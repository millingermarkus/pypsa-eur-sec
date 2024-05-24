# coding: utf-8

import pypsa
import re
import os

import pandas as pd
import numpy as np
import xarray as xr
import networkx as nx

from itertools import product
from scipy.stats import beta
from vresutils.costdata import annuity

from build_energy_totals import build_eea_co2, build_eurostat_co2, build_co2_totals
from helper import override_component_attrs, generate_periodic_profiles, update_config_with_sector_opts

from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from networkx.algorithms import complement
from pypsa.geo import haversine_pts

from networkx.algorithms.connectivity.edge_augmentation import k_edge_augmentation
from networkx.algorithms import complement
from pypsa.geo import haversine_pts

import logging

logger = logging.getLogger(__name__)

from types import SimpleNamespace
spatial = SimpleNamespace()


def define_spatial(nodes, options):
    """
    Namespace for spatial

    Parameters
    ----------
    nodes : list-like
    """

    global spatial

    spatial.nodes = nodes

    # biomass

    spatial.biomass = SimpleNamespace()

    # if options["biomass_transport"]:
    spatial.biomass.nodes = nodes + " solid biomass"
    spatial.biomass.locations = nodes
    spatial.biomass.industry = nodes + " solid biomass for industry"
    spatial.biomass.industry_cc = nodes + " solid biomass for industry CC"
    # else:
    #     spatial.biomass.nodes = nodes + " solid biomass"# = ["EU solid biomass"]
    #     spatial.biomass.locations = nodes #["EU"]
    #     spatial.biomass.industry = ["solid biomass for industry"]
    #     spatial.biomass.industry_cc = ["solid biomass for industry CC"]

    spatial.biomass.df = pd.DataFrame(vars(spatial.biomass), index=nodes)

    # co2

    spatial.co2 = SimpleNamespace()

    if options["co2_network"]:
        spatial.co2.nodes = nodes + " co2 stored"
        spatial.co2.locations = nodes
        spatial.co2.vents = nodes + " co2 vent"
    else:
        spatial.co2.nodes = ["co2 stored"]
        spatial.co2.locations = ["EU"]
        spatial.co2.vents = ["co2 vent"]

    spatial.co2.df = pd.DataFrame(vars(spatial.co2), index=nodes)

    # gas

    spatial.gas = SimpleNamespace()

    if options["gas_network"]:
        spatial.gas.nodes = nodes + " gas"
        spatial.gas.locations = nodes
        spatial.gas.biogas = nodes + " digestible biomass"
        spatial.gas.industry = nodes + " gas for industry"
        spatial.gas.industry_cc = nodes + " gas for industry CC"
        spatial.gas.biogas_to_gas = nodes + " digestible biomass to gas"
    else:
        spatial.gas.nodes = ["EU gas"]
        spatial.gas.locations = ["EU"]
        spatial.gas.biogas = ["EU biogas"]
        spatial.gas.industry = ["gas for industry"]
        spatial.gas.industry_cc = ["gas for industry CC"]
        spatial.gas.biogas_to_gas = ["EU biogas to gas"]

    spatial.gas.df = pd.DataFrame(vars(spatial.gas), index=nodes)

    # ammonia

    if options.get('ammonia'):
        spatial.ammonia = SimpleNamespace()
        if options.get("ammonia") == "regional":
            spatial.ammonia.nodes = nodes + " NH3"
            spatial.ammonia.locations = nodes
        else:
            spatial.ammonia.nodes = ["EU NH3"]
            spatial.ammonia.locations = ["EU"]

        spatial.ammonia.df = pd.DataFrame(vars(spatial.ammonia), index=nodes)

    # oil
    spatial.oil = SimpleNamespace()
    spatial.oil.nodes = ["EU oil"]
    spatial.oil.locations = ["EU"]

    # uranium
    spatial.uranium = SimpleNamespace()
    spatial.uranium.nodes = ["EU uranium"]
    spatial.uranium.locations = ["EU"]

    # coal
    spatial.coal = SimpleNamespace()
    spatial.coal.nodes = ["EU coal"]
    spatial.coal.locations = ["EU"]

    # lignite
    spatial.lignite = SimpleNamespace()
    spatial.lignite.nodes = ["EU lignite"]
    spatial.lignite.locations = ["EU"]

    return spatial

from types import SimpleNamespace
spatial = SimpleNamespace()


def emission_sectors_from_opts(opts):
    sectors = ["electricity"]
    if "T" in opts:
        sectors += [
            "rail non-elec",
            "road non-elec"
        ]
    if "H" in opts:
        sectors += [
            "residential non-elec",
            "services non-elec"
        ]
    if "I" in opts:
        sectors += [
            "industrial non-elec",
            "industrial processes",
            "domestic aviation",
            "international aviation",
            "domestic navigation",
            "international navigation"
        ]
    if "A" in opts:
        sectors += [
            "agriculture"
        ]

    return sectors


def get(item, investment_year=None):
    """Check whether item depends on investment year"""
    if isinstance(item, dict):
        return item[investment_year]
    else:
        return item


def co2_emissions_year(countries, input_eurostat, opts, emissions_scope, report_year, year):
    """
    Calculate CO2 emissions in one specific year (e.g. 1990 or 2018).
    """
    emissions_scope = snakemake.config["energy"]["emissions"]
    eea_co2 = build_eea_co2(snakemake.input.co2, year, emissions_scope)

    # TODO: read Eurostat data from year > 2014
    # this only affects the estimation of CO2 emissions for BA, RS, AL, ME, MK
    report_year = snakemake.config["energy"]["eurostat_report_year"]
    if year > 2014:
        eurostat_co2 = build_eurostat_co2(input_eurostat, countries, report_year, year=2014)
    else:
        eurostat_co2 = build_eurostat_co2(input_eurostat, countries, report_year, year)

    co2_totals = build_co2_totals(countries, eea_co2, eurostat_co2)

    sectors = emission_sectors_from_opts(opts)

    co2_emissions = co2_totals.loc[countries, sectors].sum().sum()

    # convert MtCO2 to GtCO2
    co2_emissions *= 0.

    return co2_emissions


# TODO: move to own rule with sector-opts wildcard?
def build_carbon_budget(o, input_eurostat, fn, emissions_scope, report_year):
    """
    Distribute carbon budget following beta or exponential transition path.
    """
    # opts?

    if "be" in o:
        # beta decay
        carbon_budget = float(o[o.find("cb") + 2:o.find("be")])
        be = float(o[o.find("be") + 2:])
    if "ex" in o:
        # exponential decay
        carbon_budget = float(o[o.find("cb") + 2:o.find("ex")])
        r = float(o[o.find("ex") + 2:])

    countries = n.buses.country.dropna().unique()

    e_1990 = co2_emissions_year(countries, input_eurostat, opts, emissions_scope,
                                report_year, year=1990)

    #emissions at the beginning of the path (last year available 2018)
    e_0 = co2_emissions_year(countries, input_eurostat, opts, emissions_scope,
                             report_year,year=2018)

    planning_horizons = snakemake.config['scenario']['planning_horizons']
    t_0 = planning_horizons[0]

    if "be" in o:
        # final year in the path
        t_f = t_0 + (2 * carbon_budget / e_0).round(0)

        def beta_decay(t):
            cdf_term = (t - t_0) / (t_f - t_0)
            return (e_0 / e_1990) * (1 - beta.cdf(cdf_term, be, be))

        # emissions (relative to 1990)
        co2_cap = pd.Series({t: beta_decay(t) for t in planning_horizons}, name=o)

    if "ex" in o:
        T = carbon_budget / e_0
        m = (1 + np.sqrt(1 + r * T)) / T

        def exponential_decay(t):
            return (e_0 / e_1990) * (1 + (m + r) * (t - t_0)) * np.exp(-m * (t - t_0))

        co2_cap = pd.Series({t: exponential_decay(t) for t in planning_horizons}, name=o)

    # TODO log in Snakefile
    csvs_folder = fn.rsplit("/", 1)[0]
    if not os.path.exists(csvs_folder):
        os.makedirs(csvs_folder)
    co2_cap.to_csv(fn, float_format='%.3f')


def add_lifetime_wind_solar(n, costs):
    """Add lifetime for solar and wind generators."""
    for carrier in ['solar', 'onwind', 'offwind']:
        gen_i = n.generators.index.str.contains(carrier)
        n.generators.loc[gen_i, "lifetime"] = costs.at[carrier, 'lifetime']


def haversine(p):
    coord0 = n.buses.loc[p.bus0, ['x', 'y']].values
    coord1 = n.buses.loc[p.bus1, ['x', 'y']].values
    return 1.5 * haversine_pts(coord0, coord1)


def create_network_topology(n, prefix, carriers=["DC"], connector=" -> ", bidirectional=True):
    """
    Create a network topology from transmission lines and link carrier selection.

    Parameters
    ----------
    n : pypsa.Network
    prefix : str
    carriers : list-like
    connector : str
    bidirectional : bool, default True
        True: one link for each connection
        False: one link for each connection and direction (back and forth)

    Returns
    -------
    pd.DataFrame with columns bus0, bus1, length, underwater_fraction
    """

    ln_attrs = ["bus0", "bus1", "length"]
    lk_attrs = ["bus0", "bus1", "length", "underwater_fraction"]
    lk_attrs = n.links.columns.intersection(lk_attrs)

    candidates = pd.concat([
        n.lines[ln_attrs],
        n.links.loc[n.links.carrier.isin(carriers), lk_attrs]
    ]).fillna(0)

    # base network topology purely on location not carrier
    candidates["bus0"] = candidates.bus0.map(n.buses.location)
    candidates["bus1"] = candidates.bus1.map(n.buses.location)

    positive_order = candidates.bus0 < candidates.bus1
    candidates_p = candidates[positive_order]
    swap_buses = {"bus0": "bus1", "bus1": "bus0"}
    candidates_n = candidates[~positive_order].rename(columns=swap_buses)
    candidates = pd.concat([candidates_p, candidates_n])

    def make_index(c):
        return prefix + c.bus0 + connector + c.bus1

    topo = candidates.groupby(["bus0", "bus1"], as_index=False).mean()
    topo.index = topo.apply(make_index, axis=1)

    if not bidirectional:
        topo_reverse = topo.copy()
        topo_reverse.rename(columns=swap_buses, inplace=True)
        topo_reverse.index = topo_reverse.apply(make_index, axis=1)
        topo = pd.concat([topo, topo_reverse])

    return topo


# TODO merge issue with PyPSA-Eur
def update_wind_solar_costs(n, costs):
    """
    Update costs for wind and solar generators added with pypsa-eur to those
    cost in the planning year
    """

    # NB: solar costs are also manipulated for rooftop
    # when distribution grid is inserted
    n.generators.loc[n.generators.carrier == 'solar', 'capital_cost'] = costs.at['solar-utility', 'fixed']

    n.generators.loc[n.generators.carrier == 'onwind', 'capital_cost'] = costs.at['onwind', 'fixed']

    # for offshore wind, need to calculated connection costs

    # assign clustered bus
    # map initial network -> simplified network
    busmap_s = pd.read_csv(snakemake.input.busmap_s, index_col=0).squeeze()
    busmap_s.index = busmap_s.index.astype(str)
    busmap_s = busmap_s.astype(str)
    # map simplified network -> clustered network
    busmap = pd.read_csv(snakemake.input.busmap, index_col=0).squeeze()
    busmap.index = busmap.index.astype(str)
    busmap = busmap.astype(str)
    # map initial network -> clustered network
    clustermaps = busmap_s.map(busmap)

    # code adapted from pypsa-eur/scripts/add_electricity.py
    for connection in ['dc', 'ac']:
        tech = "offwind-" + connection
        profile = snakemake.input['profile_offwind_' + connection]
        with xr.open_dataset(profile) as ds:
            underwater_fraction = ds['underwater_fraction'].to_pandas()
            connection_cost = (snakemake.config['costs']['lines']['length_factor'] *
                               ds['average_distance'].to_pandas() *
                               (underwater_fraction *
                                costs.at[tech + '-connection-submarine', 'fixed'] +
                                (1. - underwater_fraction) *
                                costs.at[tech + '-connection-underground', 'fixed']))

            # convert to aggregated clusters with weighting
            weight = ds['weight'].to_pandas()

            # e.g. clusters == 37m means that VRE generators are left
            # at clustering of simplified network, but that they are
            # connected to 37-node network
            if snakemake.wildcards.clusters[-1:] == "m":
                genmap = busmap_s
            else:
                genmap = clustermaps

            connection_cost = (connection_cost * weight).groupby(genmap).sum() / weight.groupby(genmap).sum()

            capital_cost = (costs.at['offwind', 'fixed'] +
                            costs.at[tech + '-station', 'fixed'] +
                            connection_cost)

            logger.info("Added connection cost of {:0.0f}-{:0.0f} Eur/MW/a to {}"
                        .format(connection_cost[0].min(), connection_cost[0].max(), tech))

            n.generators.loc[n.generators.carrier == tech, 'capital_cost'] = capital_cost.rename(
                index=lambda node: node + ' ' + tech)


def add_carrier_buses(n, carrier, nodes=None):
    """
    Add buses to connect e.g. coal, nuclear and oil plants
    """

    if nodes is None:
        nodes = vars(spatial)[carrier].nodes
    location = vars(spatial)[carrier].locations

    # skip if carrier already exists
    if carrier in n.carriers.index:
        return

    if not isinstance(nodes, pd.Index):
        nodes = pd.Index(nodes)

    n.add("Carrier", carrier)

    unit = "MWh_LHV" if carrier == "gas" else "MWh_th"

    n.madd("Bus",
        nodes,
        location=location,
        carrier=carrier,
        unit=unit
    )

    #capital cost could be corrected to e.g. 0.2 EUR/kWh * annuity and O&M
    n.madd("Store",
        nodes + " Store",
        bus=nodes,
        e_nom_extendable=True,
        e_cyclic=True,
        carrier=carrier,
    )

    fossils = ['coal', 'gas', 'oil', 'lignite']
    if 'nofossil' in opts and carrier in fossils:
        print('Not adding fossil ', carrier)
        extendable = False
    else:
        print('Adding fossil ', carrier)
        extendable = True

    n.madd("Generator",
        nodes,
        bus=nodes,
        p_nom_extendable=extendable,
        carrier=carrier,
        marginal_cost=costs.at[carrier, 'fuel']
    )


# TODO: PyPSA-Eur merge issue
def remove_elec_base_techs(n):
    """remove conventional generators (e.g. OCGT) and storage units (e.g. batteries and H2)
    from base electricity-only network, since they're added here differently using links
    """

    for c in n.iterate_components(snakemake.config["pypsa_eur"]):
        to_keep = snakemake.config["pypsa_eur"][c.name]
        to_remove = pd.Index(c.df.carrier.unique()).symmetric_difference(to_keep)
        print("Removing", c.list_name, "with carrier", to_remove)
        names = c.df.index[c.df.carrier.isin(to_remove)]
        n.mremove(c.name, names)
        n.carriers.drop(to_remove, inplace=True, errors="ignore")


# TODO: PyPSA-Eur merge issue
def remove_non_electric_buses(n):
    """
    remove buses from pypsa-eur with carriers which are not AC buses
    """
    print("drop buses from PyPSA-Eur with carrier: ", n.buses[~n.buses.carrier.isin(["AC", "DC"])].carrier.unique())
    n.buses = n.buses[n.buses.carrier.isin(["AC", "DC"])]


def patch_electricity_network(n):
    remove_elec_base_techs(n)
    remove_non_electric_buses(n)
    update_wind_solar_costs(n, costs)
    n.loads["carrier"] = "electricity"
    n.buses["location"] = n.buses.index
    n.buses["unit"] = "MWh_el"
    # remove trailing white space of load index until new PyPSA version after v0.18.
    n.loads.rename(lambda x: x.strip(), inplace=True)
    n.loads_t.p_set.rename(lambda x: x.strip(), axis=1, inplace=True)


def add_co2_tracking(n, options, carbon_sequestration_cost):
    # minus sign because opposite to how fossil fuels used:
    # CH4 burning puts CH4 down, atmosphere up
    n.add("Carrier", "co2",
          co2_emissions=-1.)

    # this tracks CO2 in the atmosphere
    n.add("Bus",
        "co2 atmosphere",
        location="EU",
        carrier="co2",
        unit="t_co2"
    )

    # can also be negative
    n.add("Store",
          "co2 atmosphere",
          e_nom_extendable=True,
          e_min_pu=-1,
          carrier="co2",
          bus="co2 atmosphere"
          )

    # this tracks CO2 stored, e.g. underground
    n.madd("Bus",
        spatial.co2.nodes,
        location=spatial.co2.locations,
        carrier="co2 stored",
        unit="t_co2"
    )

    print('CO2 sequestration cost: ', carbon_sequestration_cost)

    n.madd("Store",
        spatial.co2.nodes,
        e_nom_extendable=True,
        e_nom_max=np.inf,
        capital_cost=carbon_sequestration_cost,
        carrier="co2 stored",
        bus=spatial.co2.nodes
    )

    n.madd("Link",
        spatial.co2.vents,
        bus0=spatial.co2.nodes,
        bus1="co2 atmosphere",
        carrier="co2 vent",
        efficiency=1.,
        p_nom_extendable=True
    )


def add_co2_network(n, costs):

    logger.info("Adding CO2 network.")
    co2_links = create_network_topology(n, "CO2 pipeline ")

    cost_onshore = (1 - co2_links.underwater_fraction) * costs.at['CO2 pipeline', 'fixed'] * co2_links.length
    cost_submarine = co2_links.underwater_fraction * costs.at['CO2 submarine pipeline', 'fixed'] * co2_links.length
    capital_cost = cost_onshore + cost_submarine

    n.madd("Link",
        co2_links.index,
        bus0=co2_links.bus0.values + " co2 stored",
        bus1=co2_links.bus1.values + " co2 stored",
        p_min_pu=-1,
        p_nom_extendable=True,
        length=co2_links.length.values,
        capital_cost=capital_cost.values,
        carrier="CO2 pipeline",
        lifetime=costs.at['CO2 pipeline', 'lifetime']
    )

def add_KOH(n):
    #Solvent
    n.add("Carrier",
	"KOH"
	)
	    	   
    n.add("Bus",
    	"EU KOH",
	location="EU",
	carrier="KOH",
	unit="t"
	)
	    
	    
	    	   
    n.madd("Link",
	spatial.nodes + " Chlor-alkali process",
	bus0=spatial.nodes,
	bus1=spatial.nodes + " H2",
	bus2= "EU KOH",
	bus3= spatial.nodes + " >120C process steam",
	p_nom_extendable=True,
	carrier= "Chlor-alkali process",
	efficiency= 0.325, #0.2941,
		#efficiency=costs.at["electrolysis", "efficiency"],
	efficiency2 = 0.39,
	efficiency3 = - 0.86/2.88,
	capital_cost = costs.at['Chlor-alkali process', 'fixed'], #1000000,    #Euros / MW (New ranges calculated that are between 600000 and 4000000)
		#capital_cost=costs.at["electrolysis", "fixed"],
	lifetime=costs.at['Chlor-alkali process', 'lifetime'],
		#lifetime=costs.at['electrolysis', 'lifetime'],
	marginal_cost = 260/2.88 #Marginal cost, used values from DECHEMA
	)
	
    n.add("Store",
	"EU KOH store",
	bus="EU KOH",
	e_nom_extendable=True,
	e_cyclic=True,
	carrier="KOH"            #" KOH storage"    
	    	
	    	)

def add_PEI(n):
    #Solvent

    n.add("Carrier", "PEI silica sorbent")
    
    n.add("Bus",
        "EU PEI silica sorbent",
        location="EU",
        carrier="PEI silica sorbent",
        unit="t"
        )
        
    n.madd("Link",
        spatial.nodes+ " Solid sorbent production",
        bus0=spatial.nodes,
        bus1="EU PEI silica sorbent",
        bus2=spatial.nodes + " >120C process steam",     #In LCA study of the material the ragne of 35-80 degrees celsius is stated and used  (Need to assess if this is reasonable or rural/housing heat should be used instead)
        bus3="EU MEA",      #Assume initially only monoethanolamine is used and not DEA or TEA.
        bus4="EU KOH",     #Assume producing KOH and NaOH are roughly the same process and the caustic soda solution can be used for both
        efficiency= 1/0.18,   # Worst case is 0.27 kWh. Table is in kg PEI produced so converting everything to be based on MWh, best case 0.42 kWh. 0.35 is middle value that will be used
        efficiency2= - 3.8/(0.18*3.6),  #Converting 7.5 MJ to MWh and then per MWh electricity  (Took a middle value of the range)
        efficiency3= - 1/0.18, #1.42/0.35, #Best case 1.42 kg/kg PEI and worst case 2.53 kg ethanolamine per kg PEI
        efficiency4= - 1.4/0.18, #Best case 1.86 and worst case 3.32 (ALL values from LCA on solid sorbent
        p_nom_extendable=True,
        carrier="PEI silica sorbent production",
        capital_cost= costs.at['PEI production', 'fixed'], #1000000, #1018.75 from S&P estimation for plant of 320 kt/yr and 326 M$ (CONVERT TO EURO) This was for Ethylene oxide and not PEI    #1000000
        marginal_cost= costs.at['PEI production', 'VOM'], #2278/0.18, #0.70/0.18,
        lifetime=costs.at['PEI production', 'lifetime']
        )
        
        
        
    n.add("Store",
        "EU PEI silica sorbent store",
        bus="EU PEI silica sorbent",
        e_nom_extendable=True,
	e_cyclic=True,
	carrier="PEI silica sorbent"
	    	
	    	)
        


def add_dac(n, costs):
    heat_carriers = ["urban central heat"]#, "services urban decentral heat"]
    heat_buses = n.buses.index[n.buses.carrier.isin(heat_carriers)]
    locations = n.buses.location[heat_buses]

    #Assuming steam for CC to be produced by electric steam boiler
    # efficiency_el = -(costs.at['direct air capture', 'electricity-input'] + costs.at[
    #     'direct air capture', 'compression-electricity-input'] + costs.at['direct air capture', 'heat-input'] * costs.at['electric boiler steam', 'efficiency'])

    #Add steam boiler to capital cost
    # capital_cost = costs.at['direct air capture', 'fixed'] + costs.at['direct air capture', 'heat-input'] \
    #                * costs.at['electric boiler steam', 'fixed'] * costs.at['electric boiler steam', 'efficiency']
    opts = snakemake.wildcards.sector_opts.split('-')

    if "I" in opts:
        if options["industrial_steam_heat_pumps"]:
            heatbus = " <120C process steam"
        else:
            heatbus = " >120C process steam"
    else:
        heatbus = n.buses.location[heat_buses]

    efficiency_el = -(costs.at['solid direct air capture', 'electricity-input'] + costs.at['solid direct air capture', 'compression-electricity-input'])
    efficiency_th = -costs.at['solid direct air capture', 'heat-input']
    
    if options["PEI silica sorbent"]:

    	n.madd("Link",
    	    #spatial.nodes + " S-DAC",
	    heat_buses.str.replace(" heat", " DAC"),
	    bus0="co2 atmosphere",
	    bus2=spatial.co2.df.loc[locations, "nodes"].values,
	    bus1=locations.values,
	    bus3=locations.values + heatbus,
	    bus5="EU PEI silica sorbent",     #Added,  CHANGED FROM KOH TO PEI as these values are based on solid sorbent DAC
	    carrier="S-DAC",
	    capital_cost=costs.at['solid direct air capture', 'custom investment'],
	    efficiency2=min(1,options["cc_fraction"] / 0.9),#Set to a baseline of 1 for DAC when cc_fraction is >=0.9
	    efficiency=efficiency_el,
	    efficiency3=efficiency_th,
	    efficiency5= options["PEI_consumption"], #- 14e-3, #2.3e-3, #14e-3,        #Solid sorbent loss 7.5 kg / t_CO2, multiplied by 0.5 to account for the silica needed aswell. Assuming 50/50 split for silica. (Assume higher value of solvent loss is 14 kg/t_CO2
	    p_nom_extendable=True,
	    lifetime=costs.at['solid direct air capture', 'lifetime']
	       )
	  
    else: 
    
        n.madd("Link",
            #spatial.nodes + "S-DAC",
	    heat_buses.str.replace(" heat", " DAC"),
	    bus0="co2 atmosphere",
	    bus2=spatial.co2.df.loc[locations, "nodes"].values,
	    bus1=locations.values,
	    bus3=locations.values + heatbus,
	    #bus5= spatial.nodes +" PEI silica sorbent",     #Added,  CHANGED FROM KOH TO PEI as these values are based on solid sorbent DAC
	    carrier="S-DAC",
	    capital_cost=costs.at['solid direct air capture', 'custom investment'],
	    efficiency2=min(1,options["cc_fraction"] / 0.9),#Set to a baseline of 1 for DAC when cc_fraction is >=0.9
	    efficiency=efficiency_el,
	    efficiency3=efficiency_th,
	    #efficiency5= - 14e-3,        #Solid sorbent loss 7.5 kg / t_CO2, multiplied by 0.5 to account for the silica needed aswell. Assuming 50/50 split for silica. (Assume higher value of solvent loss is 14 kg/t_CO2
	    p_nom_extendable=True,
	    lifetime=costs.at['solid direct air capture', 'lifetime']
	       )
    
    if options["KOH"]:
        n.madd("Link",
            spatial.nodes + " L-DAC",
            bus0="co2 atmosphere",
            bus1=locations.values,
            bus2=spatial.co2.df.loc[locations, "nodes"].values,
            bus3=spatial.nodes + " highT industry",
            bus5="EU KOH",
            carrier = "L-DAC",
            capital_cost = costs.at['liquid direct air capture', 'custom investment'], #Should probably be 5 MEuro but discuss this. #Taken as same value as that used for Solid direct air capture (danish energy agency from 2022 found in inputs of markus technology data)
            efficiency= - (costs.at['liquid direct air capture', 'electricity-input'] + costs.at['liquid direct air capture', 'compression-electricity-input']),   #Taken the same compression of co2 as that used earlier but the electricity input from 2024 technology data from danish energy agency
            efficiency2=min(1,options["cc_fraction"] / 0.9),#Set to a baseline of 1 for DAC when cc_fraction is >=0.9,
            efficiency3=  - costs.at['liquid direct air capture', 'heat-input'],    #Heat input is based on the 2024 value for heat that was based on 2020 heat demand.
            efficiency5=  options["KOH_consumption"], #6.4e-3, #0.9e-3,    # KOH consumption based on report from (6.4 initially chosen as the higher value to be used, lower is 0.9)
            p_nom_extendable=options['KOH_p_nom_extendable'],
            lifetime=costs.at['liquid direct air capture', 'lifetime']
            )
            
    else:
         n.madd("Link",
            spatial.nodes + " L-DAC",
            bus0="co2 atmosphere",
            bus1=locations.values,
            bus2=spatial.co2.df.loc[locations, "nodes"].values,
            bus3=spatial.nodes + " highT industry",
            #bus5=spatial.nodes + " KOH",
            carrier = "L-DAC",
            capital_cost = costs.at['liquid direct air capture', 'custom investment'],  #Taken as same value as that used for Solid direct air capture (danish energy agency from 2022 found in inputs of markus technology data)
            efficiency= - (costs.at['liquid direct air capture', 'electricity-input'] + costs.at['liquid direct air capture', 'compression-electricity-input']),   #Taken the same compression of co2 as that used earlier but the electricity input from 2024 technology data from danish energy agency
            efficiency2=min(1,options["cc_fraction"] / 0.9),#Set to a baseline of 1 for DAC when cc_fraction is >=0.9,
            efficiency3= - costs.at['liquid direct air capture', 'heat-input'],    #Heat input is based on the 2024 value for heat that was based on 2020 heat demand.
            #efficiency5= - 6.4e-3,    # KOH consumption based on report from (6.4 initially chosen as the higher value to be used, lower is 0.9)
            p_nom_extendable=options['KOH_p_nom_extendable'],
            lifetime=costs.at['liquid direct air capture', 'lifetime']
            )
            
        
       
       
        



def add_co2limit(n, Nyears=1., limit=0.):

    logger.info(f"Adding CO2 budget limit as per unit of 1990 levels of {limit}")

    countries = n.buses.country.dropna().unique()

    sectors = emission_sectors_from_opts(opts)

    # convert Mt to tCO2
    co2_totals = 1e6 * pd.read_csv(snakemake.input.co2_totals_name, index_col=0)

    co2_limit = co2_totals.loc[countries, sectors].sum().sum()

    co2_limit *= limit * Nyears

    n.add("GlobalConstraint",
          "CO2Limit",
          carrier_attribute="co2_emissions",
          sense="<=",
          constant=co2_limit
          )


# TODO PyPSA-Eur merge issue

def average_every_nhours(n, offset):
    logger.info(f'Resampling the network to {offset}')
    m = n.copy(with_time=False)

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name + "_t")
        for k, df in c.pnl.items():
            if not df.empty:
                if c.list_name == "stores" and k == "e_max_pu":
                    pnl[k] = df.resample(offset).min()
                elif c.list_name == "stores" and k == "e_min_pu":
                    pnl[k] = df.resample(offset).max()
                else:
                    pnl[k] = df.resample(offset).mean()

    return m


def cycling_shift(df, steps=1):
    """Cyclic shift on index of pd.Series|pd.DataFrame by number of steps"""
    df = df.copy()
    new_index = np.roll(df.index, steps)
    df.values[:] = df.reindex(index=new_index).values
    return df


# TODO checkout PyPSA-Eur script
def prepare_costs(cost_file, USD_to_EUR, discount_rate, Nyears, lifetime):
    # set all asset costs and other parameters
    costs = pd.read_csv(cost_file, index_col=[0, 1]).sort_index()
    print(costs['value'].head(50))
    print(costs[costs['value'].isna()].head(50))
    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD"), "value"] *= USD_to_EUR

    # min_count=1 is important to generate NaNs which are then filled by fillna
    costs = costs.loc[:, "value"].unstack(level=1).groupby("technology").sum(min_count=1)
    costs = costs.fillna({"CO2 intensity": 0,
                          "FOM": 0,
                          "VOM": 0,
                          "discount rate": discount_rate,
                          "efficiency": 1,
                          "fuel": 0,
                          "investment": 0,
                          "lifetime": lifetime
                          })

    annuity_factor = lambda v: annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100
    costs["fixed"] = [annuity_factor(v) * v["investment"] * Nyears for i, v in costs.iterrows()]

    return costs


def sensitivity_costs(costs, biomass_import_price, carbon_sequestration_cost, options):

    print('Adapting costs for sensitivity analysis')
    print('Fischer-Tropsch sensitivity string: ', snakemake.wildcards.ft_s)
    print('Electrofuel sensitivity string: ', snakemake.wildcards.efu_s)
    print('Hydrogen sensitivity string: ', snakemake.wildcards.h2_s)
    print('CC sensitivity string: ', snakemake.wildcards.cc_s)
    print('CS sensitivity string: ', snakemake.wildcards.cs_s)
    print('Fossil fuel sensitivity string: ', snakemake.wildcards.fsl_s)
    print('Biomass sensitivity string: ', snakemake.wildcards.bm_s)
    print('Biomass import sensitivity string: ', snakemake.wildcards.bmim_s)

    c_in_char = 0.03
    input_CO2_intensity = costs.at['solid biomass', 'CO2 intensity']

    if 'FT0' in snakemake.wildcards.ft_s:
        # costs.at['BtL', 'efficiency'] = 0.5
        # costs.at['BtL', 'investment'] = 1500000
        costs.at['Fischer-Tropsch', 'efficiency'] = 0.9
        costs.at['Fischer-Tropsch', 'investment'] = 675000

        # costs.at['BtL', 'C in fuel'] = costs.at['BtL', 'efficiency'] * costs.at['oil', 'CO2 intensity'] / input_CO2_intensity
        # costs.at['BtL', 'C stored'] = 1 - costs.at['BtL', 'C in fuel'] - c_in_char
        # costs.at['BtL', 'CO2 stored'] = input_CO2_intensity * costs.at['BtL', 'C stored']

    elif 'FT2' in snakemake.wildcards.ft_s:
        # costs.at['BtL', 'efficiency'] = 0.35
        # costs.at['BtL', 'investment'] = 2500000
        costs.at['Fischer-Tropsch', 'efficiency'] = 0.6
        costs.at['Fischer-Tropsch', 'investment'] = 1125000

        # costs.at['BtL', 'C in fuel'] = costs.at['BtL', 'efficiency'] * costs.at[
        #     'oil', 'CO2 intensity'] / input_CO2_intensity
        # costs.at['BtL', 'C stored'] = 1 - costs.at['BtL', 'C in fuel'] - c_in_char
        # costs.at['BtL', 'CO2 stored'] = input_CO2_intensity * costs.at['BtL', 'C stored']

    elif 'FT1' in snakemake.wildcards.ft_s:
        pass

    print('BtL CO2 stored: ', costs.at['BtL', 'CO2 stored'])

    if 'EFU0' in snakemake.wildcards.efu_s:
        options["cc_fraction"] = 0.95
    elif 'EFU2' in snakemake.wildcards.efu_s:
        options["cc_fraction"] = 0.85
    elif 'EFU3' in snakemake.wildcards.efu_s:
        options["cc_fraction"] = 0.8
    elif 'EFU4' in snakemake.wildcards.efu_s:
        options["cc_fraction"] = 0.75
    elif 'EFU5' in snakemake.wildcards.efu_s:
        options["cc_fraction"] = 0.7
    elif 'EFU1' in snakemake.wildcards.efu_s:
        pass


    if 'H20' in snakemake.wildcards.h2_s:
        costs.at['electrolysis', 'efficiency'] = 0.8
        costs.at['electrolysis', 'investment'] = 150000
    elif 'H22' in snakemake.wildcards.h2_s:
        costs.at['electrolysis', 'efficiency'] = 0.7
        costs.at['electrolysis', 'investment'] = 400000
    elif 'H23' in snakemake.wildcards.h2_s:
        # costs.at['electrolysis', 'efficiency'] = 0.7
        costs.at['electrolysis', 'investment'] = 600000
    elif 'H24' in snakemake.wildcards.h2_s:
        # costs.at['electrolysis', 'efficiency'] = 0.7
        costs.at['electrolysis', 'investment'] = 800000
    elif 'H21' in snakemake.wildcards.h2_s:
        pass

    if 'CC0' in snakemake.wildcards.cc_s:
        # costs.at['biomass CHP capture', 'investment'] = 1600000
        # costs.at['cement capture', 'investment'] = 1800000
        costs.at['direct air capture', 'investment'] = 7000000
    elif 'CC2' in snakemake.wildcards.cc_s:
        # costs.at['biomass CHP capture', 'investment'] = 2400000
        # costs.at['cement capture', 'investment'] = 2600000
        costs.at['direct air capture', 'investment'] = 4000000
    elif 'CC3' in snakemake.wildcards.cc_s:
        # costs.at['biomass CHP capture', 'investment'] = 2400000
        # costs.at['cement capture', 'investment'] = 3000000
        costs.at['direct air capture', 'investment'] = 3000000
    elif 'CC4' in snakemake.wildcards.cc_s:
        # costs.at['biomass CHP capture', 'investment'] = 2700000
        # costs.at['cement capture', 'investment'] = 3000000
        costs.at['direct air capture', 'investment'] = 2000000
    elif 'CC5' in snakemake.wildcards.cc_s:
        # costs.at['biomass CHP capture', 'investment'] = 2700000
        # costs.at['cement capture', 'investment'] = 3000000
        costs.at['direct air capture', 'investment'] = 1500000
    elif 'CC1' in snakemake.wildcards.cc_s:
        pass

    if 'CS0' in snakemake.wildcards.cs_s:
        carbon_sequestration_cost = 10
    elif 'CS2' in snakemake.wildcards.cs_s:
        carbon_sequestration_cost = 50
    elif 'CS1' in snakemake.wildcards.cs_s:
        pass

    if 'FSL0' in snakemake.wildcards.fsl_s:
        costs.at["oil", 'fuel'] = 37.5
        costs.at["gas", 'fuel'] = 15
    elif 'FSL2' in snakemake.wildcards.fsl_s:
        costs.at["oil", 'fuel'] = 62.5
        costs.at["gas", 'fuel'] = 25
    elif 'FSL1' in snakemake.wildcards.fsl_s:
        pass

    if 'BM0' in snakemake.wildcards.bm_s:
        biomass_import_price = 10 * 3.6
    elif 'BM2' in snakemake.wildcards.bm_s:
        biomass_import_price = 20 * 3.6
    elif 'BM3' in snakemake.wildcards.bm_s:
        biomass_import_price = 25 * 3.6
    elif 'BM4' in snakemake.wildcards.bm_s:
        biomass_import_price = 30 * 3.6
    elif 'BM1' in snakemake.wildcards.bm_s:
        pass

    if 'IM0' in snakemake.wildcards.bmim_s:
        # biomass_import_price = 10 * 3.6
        pass
    elif 'IM2' in snakemake.wildcards.bmim_s:
        # biomass_import_price = 20 * 3.6
        pass
    elif 'IM1' in snakemake.wildcards.bmim_s:
        pass

    #Update fixed costs
    annuity_factor = lambda v: annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100
    costs["fixed"] = [annuity_factor(v) * v["investment"] * Nyears for i, v in costs.iterrows()]

    print('BtL investment: ', costs.at['BtL', 'investment'])
    print('Electrofuel investment: ', costs.at['Fischer-Tropsch', 'investment'])
    print('Biomass import price: ', biomass_import_price)

    return costs, biomass_import_price, carbon_sequestration_cost, options



def add_generation(n, costs):

    logger.info("adding electricity generation")

    nodes = pop_layout.index

    fallback = {"OCGT": "gas"}
    conventionals = options.get("conventional_generation", fallback)

    for generator, carrier in conventionals.items():

        carrier_nodes = vars(spatial)[carrier].nodes
        add_carrier_buses(n, carrier, carrier_nodes)
        n.madd("Link",
            nodes + " " + generator,
            bus0=carrier_nodes,
            bus1=nodes,
            bus3="co2 atmosphere",
            marginal_cost=costs.at[generator, 'efficiency'] * costs.at[generator, 'VOM'], #NB: VOM is per MWel
            capital_cost=costs.at[generator, 'efficiency'] * costs.at[generator, 'fixed'], #NB: fixed cost is per MWel
            p_nom_extendable=True,
            carrier=generator,
            efficiency=costs.at[generator, 'efficiency'],
            efficiency3=costs.at[carrier, 'CO2 intensity'],
            lifetime=costs.at[generator, 'lifetime']
        )


def add_ammonia(n, costs):

    logger.info("adding ammonia carrier with synthesis, cracking and storage")

    nodes = pop_layout.index

    cf_industry = snakemake.config["industry"]

    n.add("Carrier", "NH3")

    n.madd("Bus",
        spatial.ammonia.nodes,
        location=spatial.ammonia.locations,
        carrier="NH3"
    )

    #n.madd("Link",
    #    nodes,
    #    suffix=" Haber-Bosch",
    #    bus0=nodes,
    #    bus1=spatial.ammonia.nodes,
    #    bus4=nodes + " H2",
    #    p_nom_extendable=True,
    #    carrier="Haber-Bosch",
    #    efficiency=1 / (cf_industry["MWh_elec_per_tNH3_electrolysis"] / cf_industry["MWh_NH3_per_tNH3"]), # output: MW_NH3 per MW_elec
    #    efficiency4=-cf_industry["MWh_H2_per_tNH3_electrolysis"] / cf_industry["MWh_elec_per_tNH3_electrolysis"], # input: MW_H2 per MW_elec
    #    capital_cost=costs.at["Haber-Bosch", "fixed"],
    #    lifetime=costs.at["Haber-Bosch", 'lifetime']
    #)

    n.madd(
        "Link",
        nodes,
        suffix=" Haber-Bosch",
        bus0=nodes,
        bus1=spatial.ammonia.nodes,
        bus4=nodes + " H2",
        p_nom_extendable=True,
        carrier="Haber-Bosch",
        efficiency=1 /  0.2473,   #costs.at["Haber-Bosch", "electricity-input"],   #Check Technology data to see if I need to add values manually or just refer to technology data directly
        efficiency4=-1.1484 / 0.2473,   
        #costs.at["Haber-Bosch", "hydrogen-input"]
        #/ costs.at["Haber-Bosch", "electricity-input"],
        capital_cost=costs.at["Haber-Bosch", "fixed"]
        / 0.2473, #costs.at["Haber-Bosch", "electricity-input"],
        marginal_cost=costs.at["Haber-Bosch", "VOM"]
        / 0.2473, #costs.at["Haber-Bosch", "electricity-input"],
        lifetime=costs.at["Haber-Bosch", "lifetime"],
    )

    n.madd("Link",
        nodes,
        suffix=" ammonia cracker",
        bus0=spatial.ammonia.nodes,
        bus1=nodes + " H2",
        p_nom_extendable=True,
        carrier="ammonia cracker",
        efficiency=1 / cf_industry["MWh_NH3_per_MWh_H2_cracker"],  #Value is 1.46
        capital_cost=costs.at["Ammonia cracker", "fixed"] /  cf_industry["MWh_NH3_per_MWh_H2_cracker"], # given per MW_H2, also 1.46
        lifetime=costs.at['Ammonia cracker', 'lifetime']
    )

    # Ammonia Storage
    n.madd("Store",
        spatial.ammonia.nodes,
        suffix=" ammonia store",
        bus=spatial.ammonia.nodes,
        e_nom_extendable=True,
        e_cyclic=True,
        carrier="ammonia store",
        capital_cost=costs.at["NH3 (l) storage tank incl. liquefaction", "fixed"],
        lifetime=costs.at['NH3 (l) storage tank incl. liquefaction', 'lifetime']
    )
    
def add_ethylene(n,costs):
    #solvent
    n.add("Carrier", "Ethylene")
    
    n.add("Bus", 
         "EU Ethylene",
         location="EU",
         carrier="Ethylene",
         unit="t"
         )
    
    n.madd("Link",
        spatial.nodes + " Ethylene from steam cracking",  
        bus0=spatial.oil.nodes,
        bus1="EU Ethylene",    
        bus2=spatial.nodes,
        bus3="co2 atmosphere",
        p_nom_extendable=True,
        carrier="Ethylene steam cracking",
        efficiency=1/(80.6/3.6), # 80.6 is based on the fuel and also reactant part of producing ethylene unit t_hvc/MWh_naphtha using mass allocation
        efficiency2=-(0.55/3.6)/(80.6/3.6),     #Electricity scaled down to account for allocation
        
        efficiency3= 1.7/(80.6/3.6), #t_CO2/MWh_naphtha (Higher than the one proposed by Neumann in pypsa-eur)
        capital_cost = costs.at['Ethylene', 'fixed'], #1200/(80.6/3.6)*8760,      #Assume that DECHEMA value is based on 1200 Euro/t_HVC/year
        marginal_cost = 180/(80.6/3.6), #244    #Taken from levelized cost of prodution steam cracking based on average around the world. (Add source here later)
        lifetime=costs.at['Ethylene', 'lifetime']
        )
         
    n.add("Store",
    	"EU Ethylene store",
    	bus="EU Ethylene",
    	e_nom_extendable=True,
	e_cyclic=True,
    	carrier="Ethylene"
    	)

def add_EO(n):
    #Solvent   	
    n.add("Carrier", "EO")
   
    n.add("Bus",
        "EU EO",
        location="EU",
        carrier="EO",
        unit="t"
        )
       
    n.madd("Link",
        spatial.nodes + " EO prodution",
        bus0=spatial.nodes,
        bus1="EU EO",
        bus2="EU Ethylene",
        p_nom_extendable = True,
        carrier="EO production",
        efficiency=1/0.33,             #Based on Ecoinvent database paper
        efficiency2= - 0.825/0.33,
       
        capital_cost=costs.at['EO production', 'fixed'], #38237400*0.779   #Based on investment cost in report of  "Is the Liquid-Phase H2O2-Based Ethylene Oxide Process More Economical and Greener Than the Gas-Phase O2-Based Silver-Catalyzed Process?"
        lifetime=costs.at['EO production', 'lifetime'],
        marginal_cost=costs.at['EO production', 'VOM']
        ) 
       
       
    n.add("Store",
        spatial.nodes,
        suffix="EO store",
        bus="EU EO",
        e_nom_extendable=True,
        e_cyclic=True,
        carrier= "EO"
        )
       
       
       
    	
def add_MEA(n):
    #Solvent
    
    n.add("Carrier", "MEA")
    
    n.add("Bus",
        "EU MEA",
        location="EU",
        carrier= "MEA",
        unit="t"
        )
       
#    n.madd("Link",
#        spatial.nodes + " MEA production",
#        bus0=spatial.nodes + " MEA",
#        bus1=spatial.nodes + " Ethylene",
#        bus2=spatial.ammonia.nodes,                  #This might be in concentration of 20-30% instead of pure ammonia so need to check this maybe
#        bus3=spatial.nodes,
#        bus4=spatial.nodes + " >120C process steam",  #Replaces the heat required by natural gas
#        carrier="MEA",
#        p_nom_extendable=True,
#        efficiency= - 0.816,  #Allocation done by mass as suggested by ecoinvent
#        efficiency2= - 0.788,
#        efficiency3= - 0.333e-3,
#        efficiency4= - 0.5556e-3,
#        captial_cost=1000000
#        )

#    n.madd("Link",
#        spatial.nodes + " MEA production",
#        bus0=spatial.nodes,
#        bus1="EU MEA",
#        bus2=spatial.ammonia.nodes,
#        bus3="EU EO",
#        bus4=spatial.nodes+ " >120C process steam",
#        carrier="MEA productionf",
#        p_nom_extendable=True,
#        efficiency= 1/0.333,         #Based on ecoinvent paper
#        efficiency2= -0.788/0.333,
#        efficiency3= - 0.816/0.333,
#        efficiency4= - 2/3.6/0.333,
#        #capital_cost=1000000
#        marginal_cost = 382.76    #Dollar per MWh based on TAC found in destillation paper "Simulation and energy consumption evalutaion of reactive distillation process for ethanolamine production" including allocation by mass
#        )
        
        
    n.madd("Link",
        spatial.nodes + " MEA production",
        bus0=spatial.nodes+ " >120C process steam",
        bus1="EU MEA",
        bus2=spatial.ammonia.nodes,
        bus3="EU EO",
        
        carrier="MEA production",
        p_nom_extendable=True,
        efficiency= 1/1.96,         
        efficiency2= -0.25*5.166*(1/1.96), #5.166 from MWh_Nh3_per_t_Nh3 in config file
        efficiency3= - 0.75*(1/1.96),
        
        capital_cost=costs.at['MEA production', 'fixed']*(1/1.96),
        
        marginal_cost=costs.at['MEA production', 'VOM'],
        lifetime=costs.at['MEA production','lifetime'] 
        #marginal_cost = 382.76    #Dollar per MWh based on TAC found in destillation paper "Simulation and energy consumption evalutaion of reactive distillation process for ethanolamine production" including allocation by mass
        )
        

    n.add("Store",
    	"EU MEA",
    	bus="EU MEA",
    	e_nom_extendable=True,
	e_cyclic=True,
    	carrier="MEA"
    	
    	)
         


def add_wave(n, wave_cost_factor):
    # TODO: handle in Snakefile
    wave_fn = "data/WindWaveWEC_GLTB.xlsx"

    # in kW
    capacity = pd.Series({"Attenuator": 750,
                          "F2HB": 1000,
                          "MultiPA": 600})

    # in EUR/MW
    annuity_factor = annuity(25, 0.07) + 0.03
    costs = 1e6 * wave_cost_factor * annuity_factor * pd.Series({"Attenuator": 2.5,
                                                                 "F2HB": 2,
                                                                 "MultiPA": 1.5})

    sheets = pd.read_excel(wave_fn, sheet_name=["FirthForth", "Hebrides"],
                           usecols=["Attenuator", "F2HB", "MultiPA"],
                           index_col=0, skiprows=[0], parse_dates=True)

    wave = pd.concat([sheets[l].divide(capacity, axis=1) for l in locations],
                     keys=locations,
                     axis=1)

    for wave_type in costs.index:
        n.add("Generator",
              "Hebrides " + wave_type,
              bus="GB4 0",  # TODO this location is hardcoded
              p_nom_extendable=True,
              carrier="wave",
              capital_cost=costs[wave_type],
              p_max_pu=wave["Hebrides", wave_type]
              )


def insert_electricity_distribution_grid(n, costs):
    # TODO pop_layout?
    # TODO options?

    print("Inserting electricity distribution grid with investment cost factor of",
          options['electricity_distribution_grid_cost_factor'])

    nodes = pop_layout.index

    cost_factor = options['electricity_distribution_grid_cost_factor']

    n.madd("Bus",
        nodes + " low voltage",
        location=nodes,
        carrier="low voltage",
        unit="MWh_el"
    )

    n.madd("Link",
        nodes + " electricity distribution grid",
        bus0=nodes,
        bus1=nodes + " low voltage",
        p_nom_extendable=True,
        p_min_pu=-1,
        carrier="electricity distribution grid",
        efficiency=1,
        lifetime=costs.at['electricity distribution grid', 'lifetime'],
        capital_cost=costs.at['electricity distribution grid', 'fixed'] * cost_factor
    )

    # this catches regular electricity load and "industry electricity" and
    # "agriculture machinery electric" and "agriculture electricity"
    loads = n.loads.index[n.loads.carrier.str.contains("electric")]
    n.loads.loc[loads, "bus"] += " low voltage"

    bevs = n.links.index[n.links.carrier == "BEV charger"]
    n.links.loc[bevs, "bus0"] += " low voltage"

    v2gs = n.links.index[n.links.carrier == "V2G"]
    n.links.loc[v2gs, "bus1"] += " low voltage"

    hps = n.links.index[n.links.carrier.str.contains("heat pump")]
    n.links.loc[hps, "bus0"] += " low voltage"

    rh = n.links.index[n.links.carrier.str.contains("resistive heater")]
    n.links.loc[rh, "bus0"] += " low voltage"

    mchp = n.links.index[n.links.carrier.str.contains("micro gas")]
    n.links.loc[mchp, "bus1"] += " low voltage"

    # set existing solar to cost of utility cost rather the 50-50 rooftop-utility
    solar = n.generators.index[n.generators.carrier == "solar"]
    n.generators.loc[solar, "capital_cost"] = costs.at['solar-utility', 'fixed']
    if snakemake.wildcards.clusters[-1:] == "m":
        simplified_pop_layout = pd.read_csv(snakemake.input.simplified_pop_layout, index_col=0)
        pop_solar = simplified_pop_layout.total.rename(index=lambda x: x + " solar")
    else:
        pop_solar = pop_layout.total.rename(index=lambda x: x + " solar")

    # add max solar rooftop potential assuming 0.1 kW/m2 and 10 m2/person,
    # i.e. 1 kW/person (population data is in thousands of people) so we get MW
    potential = 0.1 * 10 * pop_solar

    n.madd("Generator",
        solar,
        suffix=" rooftop",
        bus=n.generators.loc[solar, "bus"] + " low voltage",
        carrier="solar rooftop",
        p_nom_extendable=True,
        p_nom_max=potential,
        marginal_cost=n.generators.loc[solar, 'marginal_cost'],
        capital_cost=costs.at['solar-rooftop', 'fixed'],
        efficiency=n.generators.loc[solar, 'efficiency'],
        p_max_pu=n.generators_t.p_max_pu[solar],
        lifetime=costs.at['solar-rooftop', 'lifetime']
    )

    n.add("Carrier", "home battery")

    n.madd("Bus",
        nodes + " home battery",
        location=nodes,
        carrier="home battery",
        unit="MWh_el"
    )

    n.madd("Store",
           nodes + " home battery",
           bus=nodes + " home battery",
           e_cyclic=True,
           e_nom_extendable=True,
           carrier="home battery",
           capital_cost=costs.at['home battery storage', 'fixed'],
           lifetime=costs.at['battery storage', 'lifetime']
           )

    n.madd("Link",
           nodes + " home battery charger",
           bus0=nodes + " low voltage",
           bus1=nodes + " home battery",
           carrier="home battery charger",
           efficiency=costs.at['battery inverter', 'efficiency'] ** 0.5,
           capital_cost=costs.at['home battery inverter', 'fixed'],
           p_nom_extendable=True,
           lifetime=costs.at['battery inverter', 'lifetime']
           )

    n.madd("Link",
           nodes + " home battery discharger",
           bus0=nodes + " home battery",
           bus1=nodes + " low voltage",
           carrier="home battery discharger",
           efficiency=costs.at['battery inverter', 'efficiency'] ** 0.5,
           marginal_cost=options['marginal_cost_storage'],
           p_nom_extendable=True,
           lifetime=costs.at['battery inverter', 'lifetime']
           )


def insert_gas_distribution_costs(n, costs):
    # TODO options?

    f_costs = options['gas_distribution_grid_cost_factor']

    print("Inserting gas distribution grid with investment cost factor of", f_costs)

    capital_cost = costs.loc['electricity distribution grid']["fixed"] * f_costs

    # gas boilers
    gas_b = n.links.index[n.links.carrier.str.contains("gas boiler") &
                          (~n.links.carrier.str.contains("urban central"))]
    n.links.loc[gas_b, "capital_cost"] += capital_cost

    # micro CHPs
    mchp = n.links.index[n.links.carrier.str.contains("micro gas")]
    n.links.loc[mchp, "capital_cost"] += capital_cost

    # TODO: Research methane grid costs and implement!
    # biogas
    biogas = n.links.index[n.links.carrier.str.contains("biogas")]
    n.links.loc[biogas, "capital_cost"] += costs.loc['electricity distribution grid']["fixed"] * f_costs

    biosng = n.links.index[n.links.carrier.str.contains("solid biomass to gas")]
    n.links.loc[biosng, "capital_cost"] += costs.loc['electricity distribution grid']["fixed"] * f_costs

    # lowT steam methane
    mchp = n.links.index[n.links.carrier.str.contains("methane for lowT industry")]
    n.links.loc[mchp, "capital_cost"] += costs.loc['electricity distribution grid']["fixed"] * f_costs


def add_electricity_grid_connection(n, costs):
    carriers = ["onwind", "solar"]

    gens = n.generators.index[n.generators.carrier.isin(carriers)]

    n.generators.loc[gens, "capital_cost"] += costs.at['electricity grid connection', 'fixed']


def add_storage_and_grids(n, costs):

    logger.info("Add hydrogen storage")

    nodes = pop_layout.index

    n.add("Carrier", "H2")

    n.madd("Bus",
        nodes + " H2",
        location=nodes,
        carrier="H2",
        unit="MWh_LHV"
    )

    n.madd("Link",
           nodes + " H2 Electrolysis",
           bus1=nodes + " H2",
           bus0=nodes,
           p_nom_extendable=True,
           carrier="H2 Electrolysis",
           efficiency=costs.at["electrolysis", "efficiency"],
           capital_cost=costs.at["electrolysis", "fixed"],
           lifetime=costs.at['electrolysis', 'lifetime']
           )

    n.madd("Link",
           nodes + " H2 Fuel Cell",
           bus0=nodes + " H2",
           bus1=nodes,
           p_nom_extendable=True,
           carrier="H2 Fuel Cell",
           efficiency=costs.at["fuel cell", "efficiency"],
           capital_cost=costs.at["fuel cell", "fixed"] * costs.at["fuel cell", "efficiency"],
           # NB: fixed cost is per MWel
           lifetime=costs.at['fuel cell', 'lifetime']
           )

    cavern_types = snakemake.config["sector"]["hydrogen_underground_storage_locations"]
    h2_caverns = pd.read_csv(snakemake.input.h2_cavern, index_col=0)

    if not h2_caverns.empty and options['hydrogen_underground_storage']:

        h2_caverns = h2_caverns[cavern_types].sum(axis=1)

        # only use sites with at least 2 TWh potential
        h2_caverns = h2_caverns[h2_caverns > 2]

        # convert TWh to MWh
        h2_caverns = h2_caverns * 1e6

        # clip at 1000 TWh for one location
        h2_caverns.clip(upper=1e9, inplace=True)

        logger.info("Add hydrogen underground storage")

        h2_capital_cost = costs.at["hydrogen storage underground", "fixed"]

        n.madd("Store",
            h2_caverns.index + " H2 Store",
            bus=h2_caverns.index + " H2",
            e_nom_extendable=True,
            e_nom_max=h2_caverns.values,
            e_cyclic=True,
            carrier="H2 Store",
            capital_cost=h2_capital_cost,
            lifetime=costs.at["hydrogen storage underground", "lifetime"]
        )

    # hydrogen stored overground (where not already underground)
    h2_capital_cost = costs.at["hydrogen storage tank type 1 including compressor", "fixed"]
    nodes_overground = h2_caverns.index.symmetric_difference(nodes)

    n.madd("Store",
           nodes_overground + " H2 Store",
           bus=nodes_overground + " H2",
           e_nom_extendable=True,
           e_cyclic=True,
           carrier="H2 Store",
           capital_cost=h2_capital_cost
           )

    if options["gas_network"] or options["H2_retrofit"]:

        fn = snakemake.input.clustered_gas_network
        gas_pipes = pd.read_csv(fn, index_col=0)

    if options["gas_network"]:

        logger.info("Add natural gas infrastructure, incl. LNG terminals, production and entry-points.")

        if options["H2_retrofit"]:
            gas_pipes["p_nom_max"] = gas_pipes.p_nom
            gas_pipes["p_nom_min"] = 0.
            # 0.1 EUR/MWkm/a to prefer decommissioning to address degeneracy
            gas_pipes["capital_cost"] = 0.1 * gas_pipes.length
        else:
            gas_pipes["p_nom_max"] = np.inf
            gas_pipes["p_nom_min"] = gas_pipes.p_nom
            gas_pipes["capital_cost"] = gas_pipes.length * costs.at['CH4 (g) pipeline', 'fixed']

        n.madd("Link",
            gas_pipes.index,
            bus0=gas_pipes.bus0 + " gas",
            bus1=gas_pipes.bus1 + " gas",
            p_min_pu=gas_pipes.p_min_pu,
            p_nom=gas_pipes.p_nom,
            p_nom_extendable=True,
            p_nom_max=gas_pipes.p_nom_max,
            p_nom_min=gas_pipes.p_nom_min,
            length=gas_pipes.length,
            capital_cost=gas_pipes.capital_cost,
            tags=gas_pipes.name,
            carrier="gas pipeline",
            lifetime=costs.at['CH4 (g) pipeline', 'lifetime']
        )

        # remove fossil generators where there is neither
        # production, LNG terminal, nor entry-point beyond system scope

        fn = snakemake.input.gas_input_nodes_simplified
        gas_input_nodes = pd.read_csv(fn, index_col=0)

        unique = gas_input_nodes.index.unique()
        gas_i = n.generators.carrier == 'gas'
        internal_i = ~n.generators.bus.map(n.buses.location).isin(unique)

        remove_i = n.generators[gas_i & internal_i].index
        n.generators.drop(remove_i, inplace=True)

        p_nom = gas_input_nodes.sum(axis=1).rename(lambda x: x + " gas")
        n.generators.loc[gas_i, "p_nom_extendable"] = False
        n.generators.loc[gas_i, "p_nom"] = p_nom

        # add candidates for new gas pipelines to achieve full connectivity

        G = nx.Graph()

        gas_buses = n.buses.loc[n.buses.carrier=='gas', 'location']
        G.add_nodes_from(np.unique(gas_buses.values))

        sel = gas_pipes.p_nom > 1500
        attrs = ["bus0", "bus1", "length"]
        G.add_weighted_edges_from(gas_pipes.loc[sel, attrs].values)

        # find all complement edges
        complement_edges = pd.DataFrame(complement(G).edges, columns=["bus0", "bus1"])
        complement_edges["length"] = complement_edges.apply(haversine, axis=1)

        # apply k_edge_augmentation weighted by length of complement edges
        k_edge = options.get("gas_network_connectivity_upgrade", 3)
        augmentation = list(k_edge_augmentation(G, k_edge, avail=complement_edges.values))

        if augmentation:

            new_gas_pipes = pd.DataFrame(augmentation, columns=["bus0", "bus1"])
            new_gas_pipes["length"] = new_gas_pipes.apply(haversine, axis=1)

            new_gas_pipes.index = new_gas_pipes.apply(
                lambda x: f"gas pipeline new {x.bus0} <-> {x.bus1}", axis=1)

            n.madd("Link",
                new_gas_pipes.index,
                bus0=new_gas_pipes.bus0 + " gas",
                bus1=new_gas_pipes.bus1 + " gas",
                p_min_pu=-1, # new gas pipes are bidirectional
                p_nom_extendable=True,
                length=new_gas_pipes.length,
                capital_cost=new_gas_pipes.length * costs.at['CH4 (g) pipeline', 'fixed'],
                carrier="gas pipeline new",
                lifetime=costs.at['CH4 (g) pipeline', 'lifetime']
            )

    if options["H2_retrofit"]:

        logger.info("Add retrofitting options of existing CH4 pipes to H2 pipes.")

        fr = "gas pipeline"
        to = "H2 pipeline retrofitted"
        h2_pipes = gas_pipes.rename(index=lambda x: x.replace(fr, to))

        n.madd("Link",
            h2_pipes.index,
            bus0=h2_pipes.bus0 + " H2",
            bus1=h2_pipes.bus1 + " H2",
            p_min_pu=-1.,  # allow that all H2 retrofit pipelines can be used in both directions
            p_nom_max=h2_pipes.p_nom * options["H2_retrofit_capacity_per_CH4"],
            p_nom_extendable=True,
            length=h2_pipes.length,
            capital_cost=costs.at['H2 (g) pipeline repurposed', 'fixed'] * h2_pipes.length,
            tags=h2_pipes.name,
            carrier="H2 pipeline retrofitted",
            lifetime=costs.at['H2 (g) pipeline repurposed', 'lifetime']
        )

    if options.get("H2_network", True):

        logger.info("Add options for new hydrogen pipelines.")

        h2_pipes = create_network_topology(n, "H2 pipeline ", carriers=["DC", "gas pipeline"])

        # TODO Add efficiency losses
        n.madd("Link",
            h2_pipes.index,
            bus0=h2_pipes.bus0.values + " H2",
            bus1=h2_pipes.bus1.values + " H2",
            p_min_pu=-1,
            p_nom_extendable=True,
            length=h2_pipes.length.values,
            capital_cost=costs.at['H2 (g) pipeline', 'fixed'] * h2_pipes.length.values,
            carrier="H2 pipeline",
            lifetime=costs.at['H2 (g) pipeline', 'lifetime']
        )

    n.add("Carrier", "battery")

    n.madd("Bus",
        nodes + " battery",
        location=nodes,
        carrier="battery",
        unit="MWh_el"
    )

    n.madd("Store",
           nodes + " battery",
           bus=nodes + " battery",
           e_cyclic=True,
           e_nom_extendable=True,
           carrier="battery",
           capital_cost=costs.at['battery storage', 'fixed'],
           lifetime=costs.at['battery storage', 'lifetime']
           )

    n.madd("Link",
           nodes + " battery charger",
           bus0=nodes,
           bus1=nodes + " battery",
           carrier="battery charger",
           efficiency=costs.at['battery inverter', 'efficiency'] ** 0.5,
           capital_cost=costs.at['battery inverter', 'fixed'],
           p_nom_extendable=True,
           lifetime=costs.at['battery inverter', 'lifetime']
           )

    n.madd("Link",
           nodes + " battery discharger",
           bus0=nodes + " battery",
           bus1=nodes,
           carrier="battery discharger",
           efficiency=costs.at['battery inverter', 'efficiency'] ** 0.5,
           marginal_cost=options['marginal_cost_storage'],
           p_nom_extendable=True,
           lifetime=costs.at['battery inverter', 'lifetime']
           )

    if options['methanation']:
        n.madd("Link",
            spatial.nodes,
            suffix=" Sabatier",
            bus0=nodes + " H2",
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            p_nom_extendable=True,
            carrier="Sabatier",
            efficiency=costs.at["methanation", "efficiency"],
            efficiency2=-costs.at["methanation", "efficiency"] * costs.at['gas', 'CO2 intensity'],
            capital_cost=costs.at["methanation", "fixed"] * costs.at["methanation", "efficiency"],  # costs given per kW_gas
            lifetime=costs.at['methanation', 'lifetime']
        )

    if options['helmeth']:
        n.madd("Link",
            spatial.nodes,
            suffix=" helmeth",
            bus0=nodes,
            bus1=spatial.gas.nodes,
            bus2=spatial.co2.nodes,
            carrier="helmeth",
            p_nom_extendable=True,
            efficiency=costs.at["helmeth", "efficiency"],
            efficiency2=-costs.at["helmeth", "efficiency"] * costs.at['gas', 'CO2 intensity'],
            capital_cost=costs.at["helmeth", "fixed"],
            lifetime=costs.at['helmeth', 'lifetime']
        )
        
    

    if options.get('coal_cc'):

        n.madd("Link",
            spatial.nodes,
            suffix=" coal CC",
            bus0=spatial.coal.nodes,
            bus1=spatial.nodes,
            bus3="co2 atmosphere",
            bus2="co2 stored",
            bus4="EU MEA" ,   #Adding MEA consumption
            marginal_cost=costs.at['coal', 'efficiency'] * costs.at['coal', 'VOM'], #NB: VOM is per MWel
            capital_cost=costs.at['coal', 'efficiency'] * costs.at['coal', 'fixed'] + costs.at['biomass CHP capture', 'fixed'] * costs.at['coal', 'CO2 intensity'], #NB: fixed cost is per MWel
            p_nom_extendable=True,
            carrier="coal",
            efficiency=costs.at['coal', 'efficiency'],
            efficiency3=costs.at['coal', 'CO2 intensity'] * (1 - options["cc_fraction"]),
            efficiency2=costs.at['coal', 'CO2 intensity'] * options["cc_fraction"],
            efficiency4= - costs.at['coal', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"],   #MEA consumption in t per MWh_coal  (MEA higher value set to 3.98 kg /t_CO2
            lifetime=costs.at['coal','lifetime']
        )

    if options['SMR']:
        n.madd("Link",
            spatial.nodes,
            suffix=" SMR CC",
            bus0=spatial.gas.nodes,
            bus1=nodes + " H2",
            bus3="co2 atmosphere",
            bus2=spatial.co2.nodes,
            #bus5="EU MEA",     #Adding MEA consumption
            p_nom_extendable=True,
            carrier="SMR CC",
            efficiency=costs.at["SMR CC", "efficiency"],
            efficiency3=costs.at['gas', 'CO2 intensity'] * (1 - options["cc_fraction"]),
            efficiency2=costs.at['gas', 'CO2 intensity'] * options["cc_fraction"],
            #efficiency5=-options["MEA_consumption"] * costs.at["gas", "CO2 intensity"] * options["cc_fraction"], #3.98e-3,
            capital_cost=costs.at["SMR CC", "fixed"],
            lifetime=costs.at['SMR CC', 'lifetime']
        )

        n.madd("Link",
            nodes + " SMR",
            bus0=spatial.gas.nodes,
            bus1=nodes + " H2",
            bus3="co2 atmosphere",
            p_nom_extendable=True,
            carrier="SMR",
            efficiency=costs.at["SMR", "efficiency"],
            efficiency3=costs.at['gas', 'CO2 intensity'],
            capital_cost=costs.at["SMR", "fixed"],
            lifetime=costs.at['SMR', 'lifetime']
        )

def add_land_transport(n, costs):
    # TODO options?

    logger.info("Add land transport")

    transport = pd.read_csv(snakemake.input.transport_demand, index_col=0, parse_dates=True)

    # Add transport demand factor depending on the year
    transport = transport * get(options["land_transport_demand"], investment_year)

    number_cars = pd.read_csv(snakemake.input.transport_data, index_col=0)["number cars"]
    avail_profile = pd.read_csv(snakemake.input.avail_profile, index_col=0, parse_dates=True)
    dsm_profile = pd.read_csv(snakemake.input.dsm_profile, index_col=0, parse_dates=True)

    fuel_cell_share = get(options["land_transport_fuel_cell_share"], investment_year)
    electric_share = get(options["land_transport_electric_share"], investment_year)
    ice_share = 1 - fuel_cell_share - electric_share

    print("FCEV share", fuel_cell_share)
    print("EV share", electric_share)
    print("ICEV share", ice_share)

    assert ice_share >= 0, "Error, more FCEV and EV share than 1."

    nodes = pop_layout.index

    if electric_share > 0:

        # n.add("Carrier", "Li ion")

        n.madd("Bus",
            nodes,
            location=nodes,
            suffix=" EV battery",
            carrier="Li ion",
            unit="MWh_el"
        )

        p_set = electric_share * (transport[nodes] + cycling_shift(transport[nodes], 1) + cycling_shift(transport[nodes], 2)) / 3

        n.madd("Load",
            nodes,
            suffix=" land transport EV",
            bus=nodes + " EV battery",
            carrier="land transport EV",
            p_set=p_set
        )

        p_nom = number_cars * options.get("bev_charge_rate", 0.011) * electric_share

        n.madd("Link",
            nodes,
            suffix= " BEV charger",
            bus0=nodes,
            bus1=nodes + " EV battery",
            p_nom=p_nom,
            carrier="BEV charger",
            p_max_pu=avail_profile[nodes],
            efficiency=options.get("bev_charge_efficiency", 0.9),
            #These were set non-zero to find LU infeasibility when availability = 0.25
            #p_nom_extendable=True,
            #p_nom_min=p_nom,
            #capital_cost=1e6,  #i.e. so high it only gets built where necessary
        )

    if electric_share > 0 and options["v2g"]:

        n.madd("Link",
            nodes,
            suffix=" V2G",
            bus1=nodes,
            bus0=nodes + " EV battery",
            p_nom=p_nom,
            carrier="V2G",
            p_max_pu=avail_profile[nodes],
            efficiency=options.get("bev_charge_efficiency", 0.9),
        )

    if electric_share > 0 and options["bev_dsm"]:

        e_nom = number_cars * options.get("bev_energy", 0.05) * options["bev_availability"] * electric_share

        n.madd("Store",
            nodes,
            suffix=" battery storage",
            bus=nodes + " EV battery",
            carrier="battery storage",
            e_cyclic=True,
            e_nom=e_nom,
            e_max_pu=1,
            e_min_pu=dsm_profile[nodes]
        )

    if fuel_cell_share > 0:

        n.madd("Load",
            nodes,
            suffix=" land transport fuel cell",
            bus=nodes + " H2",
            carrier="land transport fuel cell",
            p_set=fuel_cell_share / options['transport_fuel_cell_efficiency'] * transport[nodes]
        )

    if ice_share > 0:

        if "oil" not in n.buses.carrier.unique():
            n.madd("Bus",
                spatial.oil.nodes,
                location=spatial.oil.locations,
                carrier="oil",
                unit="MWh_LHV"
            )

        ice_efficiency = options['transport_internal_combustion_efficiency']

        n.madd("Load",
            nodes,
            suffix=" land transport oil",
            bus=spatial.oil.nodes,
            carrier="land transport oil",
            p_set=ice_share / ice_efficiency * transport[nodes]
        )

        co2 = ice_share / ice_efficiency * transport[nodes].sum().sum() / 8760 * costs.at["oil", 'CO2 intensity']

        n.add("Load",
            "land transport oil emissions",
            bus="co2 atmosphere",
            carrier="land transport oil emissions",
            p_set=-co2
        )


def build_heat_demand(n):

    # copy forward the daily average heat demand into each hour, so it can be multipled by the intraday profile
    daily_space_heat_demand = xr.open_dataarray(snakemake.input.heat_demand_total).to_pandas().reindex(index=n.snapshots, method="ffill")

    intraday_profiles = pd.read_csv(snakemake.input.heat_profile, index_col=0)

    sectors = ["residential", "services"]
    uses = ["water", "space"]

    heat_demand = {}
    electric_heat_supply = {}
    for sector, use in product(sectors, uses):
        weekday = list(intraday_profiles[f"{sector} {use} weekday"])
        weekend = list(intraday_profiles[f"{sector} {use} weekend"])
        weekly_profile = weekday * 5 + weekend * 2
        intraday_year_profile = generate_periodic_profiles(
            daily_space_heat_demand.index.tz_localize("UTC"),
            nodes=daily_space_heat_demand.columns,
            weekly_profile=weekly_profile
        )

        if use == "space":
            heat_demand_shape = daily_space_heat_demand * intraday_year_profile
        else:
            heat_demand_shape = intraday_year_profile

        heat_demand[f"{sector} {use}"] = (heat_demand_shape/heat_demand_shape.sum()).multiply(pop_weighted_energy_totals[f"total {sector} {use}"]) * 1e6
        electric_heat_supply[f"{sector} {use}"] = (heat_demand_shape/heat_demand_shape.sum()).multiply(pop_weighted_energy_totals[f"electricity {sector} {use}"]) * 1e6

    heat_demand = pd.concat(heat_demand, axis=1)
    electric_heat_supply = pd.concat(electric_heat_supply, axis=1)

    # subtract from electricity load since heat demand already in heat_demand
    electric_nodes = n.loads.index[n.loads.carrier == "electricity"]
    n.loads_t.p_set[electric_nodes] = n.loads_t.p_set[electric_nodes] - electric_heat_supply.groupby(level=1, axis=1).sum()[electric_nodes]

    return heat_demand


def add_heat(n, costs):

    logger.info("Add heat sector")

    sectors = ["residential", "services"]


    heat_demand = build_heat_demand(n)

    nodes, dist_fraction, urban_fraction = create_nodes_for_heat_sector()

    #NB: must add costs of central heating afterwards (EUR 400 / kWpeak, 50a, 1% FOM from Fraunhofer ISE)

    # exogenously reduce space heat demand
    if options["reduce_space_heat_exogenously"]:
        dE = get(options["reduce_space_heat_exogenously_factor"], investment_year)
        print(f"assumed space heat reduction of {dE * 100} %")
        for sector in sectors:
            heat_demand[sector + " space"] = (1 - dE) * heat_demand[sector + " space"]

    heat_systems = [
        "residential rural",
        "services rural",
        "residential urban decentral",
        "services urban decentral",
        "urban central"
    ]

    cop = {
        "air": xr.open_dataarray(snakemake.input.cop_air_total).to_pandas().reindex(index=n.snapshots),
        "ground": xr.open_dataarray(snakemake.input.cop_soil_total).to_pandas().reindex(index=n.snapshots)
    }

    solar_thermal = xr.open_dataarray(snakemake.input.solar_thermal_total).to_pandas().reindex(index=n.snapshots)
    # 1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
    solar_thermal = options['solar_cf_correction'] * solar_thermal / 1e3

    for name in heat_systems:

        name_type = "central" if name == "urban central" else "decentral"

        n.add("Carrier", name + " heat")

        n.madd("Bus",
            nodes[name] + f" {name} heat",
            location=nodes[name],
            carrier=name + " heat",
            unit="MWh_th"
        )

        ## Add heat load

        for sector in sectors:
            # heat demand weighting
            if "rural" in name:
                factor = 1 - urban_fraction[nodes[name]]
            elif "urban central" in name:
                factor = dist_fraction[nodes[name]]
            elif "urban decentral" in name:
                factor = urban_fraction[nodes[name]] - \
                    dist_fraction[nodes[name]]
            else:
                raise NotImplementedError(f" {name} not in " f"heat systems: {heat_systems}")

            if sector in name:
                heat_load = heat_demand[[sector + " water", sector + " space"]].groupby(level=1, axis=1).sum()[
                    nodes[name]].multiply(factor)

        if name == "urban central":
            heat_load = heat_demand.groupby(level=1,axis=1).sum()[nodes[name]].multiply(factor * (1 + options['district_heating']['district_heating_loss']))

        n.madd("Load",
               nodes[name],
               suffix=f" {name} heat",
               bus=nodes[name] + f" {name} heat",
               carrier=name + " heat",
               p_set=heat_load
               )

        ## Add heat pumps

        heat_pump_type = "air" if "urban" in name else "ground"

        costs_name = f"{name_type} {heat_pump_type}-sourced heat pump"
        efficiency = cop[heat_pump_type][nodes[name]] if options["time_dep_hp_cop"] else costs.at[costs_name, 'efficiency']

        n.madd("Link",
               nodes[name],
               suffix=f" {name} {heat_pump_type} heat pump",
               bus0=nodes[name],
               bus1=nodes[name] + f" {name} heat",
               carrier=f"{name} {heat_pump_type} heat pump",
               efficiency=efficiency,
               capital_cost=costs.at[costs_name, 'efficiency'] * costs.at[costs_name, 'fixed'],
               p_nom_extendable=True,
               lifetime=costs.at[costs_name, 'lifetime']
               )

        if options["tes"]:

            n.add("Carrier", name + " water tanks")

            n.madd("Bus",
                nodes[name] + f" {name} water tanks",
                location=nodes[name],
                carrier=name + " water tanks",
                unit="MWh_th"
            )

            n.madd("Link",
                   nodes[name] + f" {name} water tanks charger",
                   bus0=nodes[name] + f" {name} heat",
                   bus1=nodes[name] + f" {name} water tanks",
                   efficiency=costs.at['water tank charger', 'efficiency'],
                   carrier=name + " water tanks charger",
                   p_nom_extendable=True
                   )

            n.madd("Link",
                nodes[name] + f" {name} water tanks discharger",
                bus0=nodes[name] + f" {name} water tanks",
                bus1=nodes[name] + f" {name} heat",
                carrier=name + " water tanks discharger",
                efficiency=costs.at['water tank discharger', 'efficiency'],
                p_nom_extendable=True
            )

            if isinstance(options["tes_tau"], dict):
                tes_time_constant_days = options["tes_tau"][name_type]
            else:
                logger.warning("Deprecated: a future version will require you to specify 'tes_tau' ",
                               "for 'decentral' and 'central' separately.")
                tes_time_constant_days = options["tes_tau"] if name_type == "decentral" else 180.

            n.madd("Store",
                nodes[name] + f" {name} water tanks",
                bus=nodes[name] + f" {name} water tanks",
                e_cyclic=True,
                e_nom_extendable=True,
                carrier=name + " water tanks",
                standing_loss=1 - np.exp(- 1 / 24 / tes_time_constant_days),
                capital_cost=costs.at[name_type + ' water tank storage', 'fixed'],
                lifetime=costs.at[name_type + ' water tank storage', 'lifetime']
            )

        if options["boilers"]:
            key = f"{name_type} resistive heater"

            n.madd("Link",
                   nodes[name] + f" {name} resistive heater",
                   bus0=nodes[name],
                   bus1=nodes[name] + f" {name} heat",
                   carrier=name + " resistive heater",
                   efficiency=costs.at[key, 'efficiency'],
                   capital_cost=costs.at[key, 'efficiency'] * costs.at[key, 'fixed'],
                   p_nom_extendable=True,
                   lifetime=costs.at[key, 'lifetime']
                   )

            key = f"{name_type} gas boiler"

            n.madd("Link",
                nodes[name] + f" {name} gas boiler",
                p_nom_extendable=True,
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name] + f" {name} heat",
                bus3="co2 atmosphere",
                carrier=name + " gas boiler",
                efficiency=costs.at[key, 'efficiency'],
                efficiency3=costs.at['gas', 'CO2 intensity'],
                capital_cost=costs.at[key, 'efficiency'] * costs.at[key, 'fixed'],
                lifetime=costs.at[key, 'lifetime']
            )


        if options["solar_thermal"]:
            n.add("Carrier", name + " solar thermal")

            n.madd("Generator",
                   nodes[name],
                   suffix=f" {name} solar thermal collector",
                   bus=nodes[name] + f" {name} heat",
                   carrier=name + " solar thermal",
                   p_nom_extendable=True,
                   capital_cost=costs.at[name_type + ' solar thermal', 'fixed'],
                   p_max_pu=solar_thermal[nodes[name]],
                   lifetime=costs.at[name_type + ' solar thermal', 'lifetime']
                   )

        if options["chp"] and name == "urban central":
            # add gas CHP; biomass CHP is added in biomass section
            n.madd("Link",
                nodes[name] + " urban central gas CHP",
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus4=nodes[name] + " urban central heat",
                bus3="co2 atmosphere",
                carrier="urban central gas CHP",
                p_nom_extendable=True,
                capital_cost=costs.at['central gas CHP', 'fixed'] * costs.at['central gas CHP', 'efficiency'],
                marginal_cost=costs.at['central gas CHP', 'VOM'],
                efficiency=costs.at['central gas CHP', 'efficiency'],
                efficiency4=costs.at['central gas CHP', 'efficiency'] / costs.at['central gas CHP', 'c_b'],
                efficiency3=costs.at['gas', 'CO2 intensity'],
                lifetime=costs.at['central gas CHP', 'lifetime']
            )

            n.madd("Link",
                nodes[name] + " urban central gas CHP CC",
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus4=nodes[name] + " urban central heat",
                bus3="co2 atmosphere",
                bus2=spatial.co2.df.loc[nodes[name], "nodes"].values,
                bus5 ="EU MEA",
                carrier="urban central gas CHP CC",
                p_nom_extendable=True,
                capital_cost=costs.at['central gas CHP', 'fixed']*costs.at['central gas CHP', 'efficiency'] + costs.at['biomass CHP capture', 'fixed']*costs.at['gas', 'CO2 intensity'],
                marginal_cost=costs.at['central gas CHP', 'VOM'],
                efficiency=costs.at['central gas CHP', 'efficiency'] - costs.at['gas', 'CO2 intensity'] * (costs.at['biomass CHP capture', 'electricity-input'] + costs.at['biomass CHP capture', 'compression-electricity-input']),
                efficiency4=costs.at['central gas CHP', 'efficiency'] / costs.at['central gas CHP', 'c_b'] + costs.at['gas', 'CO2 intensity'] * (costs.at['biomass CHP capture', 'heat-output'] + costs.at['biomass CHP capture', 'compression-heat-output'] - costs.at['biomass CHP capture', 'heat-input']),
                efficiency3=costs.at['gas', 'CO2 intensity'] * (1-options["cc_fraction"]),
                efficiency2=costs.at['gas', 'CO2 intensity'] * options["cc_fraction"],
                efficiency5= -costs.at['gas', 'CO2 intensity' ] * options["cc_fraction"] * options['MEA_consumption'],
                lifetime=costs.at['central gas CHP', 'lifetime']
            )

        if options["chp"] and options["hydrogen_chp"] and name == "urban central":
            n.madd("Link",
                nodes[name] + " urban central hydrogen CHP",
                bus0=nodes[name] + " H2",
                bus1=nodes[name],
                bus4=nodes[name] + " urban central heat",
                carrier="urban central hydrogen CHP",
                p_nom_extendable=True,
                capital_cost=costs.at['central hydrogen CHP', 'fixed'] * costs.at['central hydrogen CHP', 'efficiency'],
                # marginal_cost=costs.at['central hydrogen CHP', 'VOM'],
                efficiency=costs.at['central hydrogen CHP', 'efficiency'],
                efficiency4=costs.at['central hydrogen CHP', 'efficiency'] / costs.at['central hydrogen CHP', 'c_b'],
                lifetime=costs.at['central hydrogen CHP', 'lifetime']
            )

        if options["chp"] and options["micro_chp"] and name != "urban central":
            n.madd("Link",
                nodes[name] + f" {name} micro gas CHP",
                p_nom_extendable=True,
                bus0=spatial.gas.df.loc[nodes[name], "nodes"].values,
                bus1=nodes[name],
                bus4=nodes[name] + f" {name} heat",
                bus3="co2 atmosphere",
                carrier=name + " micro gas CHP",
                efficiency=costs.at['micro CHP', 'efficiency'],
                efficiency4=costs.at['micro CHP', 'efficiency-heat'],
                efficiency3=costs.at['gas', 'CO2 intensity'],
                capital_cost=costs.at['micro CHP', 'fixed'],
                lifetime=costs.at['micro CHP', 'lifetime']
            )

    if options['retrofitting']['retro_endogen']:

        logger.info("Add retrofitting endogenously")

        # resample heat demand temporal 'heat_demand_r' depending on in config
        # specified temporal resolution, to not overestimate retrofitting
        hours = list(filter(re.compile(r'^\d+h$', re.IGNORECASE).search, opts))
        if len(hours) == 0:
            hours = [n.snapshots[1] - n.snapshots[0]]
        heat_demand_r = heat_demand.resample(hours[0]).mean()

        # retrofitting data 'retro_data' with 'costs' [EUR/m^2] and heat
        # demand 'dE' [per unit of original heat demand] for each country and
        # different retrofitting strengths [additional insulation thickness in m]
        retro_data = pd.read_csv(snakemake.input.retro_cost_energy,
                                 index_col=[0, 1], skipinitialspace=True,
                                 header=[0, 1])
        # heated floor area [10^6 * m^2] per country
        floor_area = pd.read_csv(snakemake.input.floor_area, index_col=[0, 1])

        n.add("Carrier", "retrofitting")

        # share of space heat demand 'w_space' of total heat demand
        w_space = {}
        for sector in sectors:
            w_space[sector] = heat_demand_r[sector + " space"] / \
                              (heat_demand_r[sector + " space"] + heat_demand_r[sector + " water"])
        w_space["tot"] = ((heat_demand_r["services space"] +
                           heat_demand_r["residential space"]) /
                          heat_demand_r.groupby(level=[1], axis=1).sum())

        for name in n.loads[n.loads.carrier.isin([x + " heat" for x in heat_systems])].index:

            node = n.buses.loc[name, "location"]
            ct = pop_layout.loc[node, "ct"]

            # weighting 'f' depending on the size of the population at the node
            f = urban_fraction[node] if "urban" in name else (1 - urban_fraction[node])
            if f == 0:
                continue
            # get sector name ("residential"/"services"/or both "tot" for urban central)
            sec = [x if x in name else "tot" for x in sectors][0]

            # get floor aread at node and region (urban/rural) in m^2
            floor_area_node = ((pop_layout.loc[node].fraction
                                * floor_area.loc[ct, "value"] * 10 ** 6).loc[sec] * f)
            # total heat demand at node [MWh]
            demand = (n.loads_t.p_set[name].resample(hours[0])
                      .mean())

            # space heat demand at node [MWh]
            space_heat_demand = demand * w_space[sec][node]
            # normed time profile of space heat demand 'space_pu' (values between 0-1),
            # p_max_pu/p_min_pu of retrofitting generators
            space_pu = (space_heat_demand / space_heat_demand.max()).to_frame(name=node)

            # minimum heat demand 'dE' after retrofitting in units of original heat demand (values between 0-1)
            dE = retro_data.loc[(ct, sec), ("dE")]
            # get addtional energy savings 'dE_diff' between the different retrofitting strengths/generators at one node
            dE_diff = abs(dE.diff()).fillna(1 - dE.iloc[0])
            # convert costs Euro/m^2 -> Euro/MWh
            capital_cost = retro_data.loc[(ct, sec), ("cost")] * floor_area_node / \
                           ((1 - dE) * space_heat_demand.max())
            # number of possible retrofitting measures 'strengths' (set in list at config.yaml 'l_strength')
            # given in additional insulation thickness [m]
            # for each measure, a retrofitting generator is added at the node
            strengths = retro_data.columns.levels[1]

            # check that ambitious retrofitting has higher costs per MWh than moderate retrofitting
            if (capital_cost.diff() < 0).sum():
                logger.warning(f"Costs are not linear for {ct} {sec}")
                s = capital_cost[(capital_cost.diff() < 0)].index
                strengths = strengths.drop(s)

            # reindex normed time profile of space heat demand back to hourly resolution
            space_pu = space_pu.reindex(index=heat_demand.index).fillna(method="ffill")

            # add for each retrofitting strength a generator with heat generation profile following the profile of the heat demand
            for strength in strengths:
                n.madd('Generator',
                       [node],
                       suffix=' retrofitting ' + strength + " " + name[6::],
                       bus=name,
                       carrier="retrofitting",
                       p_nom_extendable=True,
                       p_nom_max=dE_diff[strength] * space_heat_demand.max(),
                       # maximum energy savings for this renovation strength
                       p_max_pu=space_pu,
                       p_min_pu=space_pu,
                       country=ct,
                       capital_cost=capital_cost[strength] * options['retrofitting']['cost_factor']
                       )


def create_nodes_for_heat_sector():
    # TODO pop_layout

    # rural are areas with low heating density and individual heating
    # urban are areas with high heating density
    # urban can be split into district heating (central) and individual heating (decentral)

    ct_urban = pop_layout.urban.groupby(pop_layout.ct).sum()
    # distribution of urban population within a country
    pop_layout["urban_ct_fraction"] = pop_layout.urban / pop_layout.ct.map(ct_urban.get)

    sectors = ["residential", "services"]

    nodes = {}
    urban_fraction = pop_layout.urban / pop_layout[["rural", "urban"]].sum(axis=1)

    for sector in sectors:
        nodes[sector + " rural"] = pop_layout.index
        nodes[sector + " urban decentral"] = pop_layout.index

    district_heat_share = pop_weighted_energy_totals["district heat share"]

    # maximum potential of urban demand covered by district heating
    central_fraction = options['district_heating']["potential"]
    # district heating share at each node
    dist_fraction_node = district_heat_share * pop_layout["urban_ct_fraction"] / pop_layout["fraction"]
    nodes["urban central"] = dist_fraction_node.index
    # if district heating share larger than urban fraction -> set urban
    # fraction to district heating share
    urban_fraction = pd.concat([urban_fraction, dist_fraction_node],
                               axis=1).max(axis=1)
    # difference of max potential and today's share of district heating
    diff = (urban_fraction * central_fraction) - dist_fraction_node
    progress = get(options["district_heating"]["progress"], investment_year)
    dist_fraction_node += diff * progress
    print(
        "The current district heating share compared to the maximum",
        f"possible is increased by a progress factor of\n{progress}",
        f"resulting in a district heating share of\n{dist_fraction_node}"
    )

    return nodes, dist_fraction_node, urban_fraction


def add_biomass(n, costs, beccs, biomass_import_price):

    logger.info("Add biomass")

    nodes = pop_layout.index

    # biomass distributed at country level - i.e. transport within country allowed
    # cts = pop_layout.ct.value_counts().index

    biomass_pot_node = pd.read_csv(snakemake.input.biomass_potentials, index_col=0)

    # need to aggregate potentials if gas not nodally resolved
    # if options["gas_network"]:
    #     biogas_potentials_spatial = biomass_potentials["biogas"].rename(index=lambda x: x + " biogas")
    # else:
    #     biogas_potentials_spatial = biomass_potentials["biogas"].sum()

    # if options["biomass_transport"]:
    #     solid_biomass_potentials_spatial = biomass_potentials["solid biomass"].rename(index=lambda x: x + " solid biomass")
    # else:
    #     solid_biomass_potentials_spatial = biomass_potentials["solid biomass"].sum()

    # potential per node distributed within country by population
    # biomass_pot_node = (biomass_potentials.loc[pop_layout.ct]
    #                     .set_index(pop_layout.index)
    #                     .mul(pop_layout.fraction, axis="index"))

    biomass_costs = pd.read_csv('resources/biomass_country_costs.csv', index_col=0)
    biomass_costs_node = (biomass_costs.loc[pop_layout.ct].set_index(pop_layout.index))
    # biomass_costs_node = biomass_costs
    # print(biomass_costs_node)
    # print(biomass_pot_node)

    # Digestible biomass
    n.add("Carrier", "digestible biomass")

    n.madd("Bus",
           nodes + " digestible biomass",
           location=nodes,
           carrier="digestible biomass")

    n.madd("Store",
           nodes + " digestible biomass",
           bus=nodes + " digestible biomass",
           carrier="digestible biomass",
           e_cyclic=True)

    # Municipal solid waste
    n.add("Carrier", "municipal solid waste")

    n.madd("Bus",
           nodes + " municipal solid waste",
           location=nodes,
           carrier="municipal solid waste")

    n.madd("Store",
           nodes + " municipal solid waste",
           bus=nodes + " municipal solid waste",
           carrier="municipal solid waste",
           e_cyclic=True)

    # Solid biomass
    n.add("Carrier", "solid biomass")

    n.madd("Bus",
           nodes + " solid biomass",
           location=nodes,
           carrier="solid biomass")

    n.madd("Store",
           nodes + " solid biomass",
           bus=nodes + " solid biomass",
           carrier="solid biomass",
           e_cyclic=True)


    biomass_types =  ["manureslurry", "sewage sludge", "straw", "forest residues", "industry wood residues", "landscape care", "municipal solid waste"]
    biomass_potential = {}
    biomass_costs = {}
    for name in biomass_types:
        biomass_potential[name] = biomass_pot_node[name].values
        biomass_costs[name] = ((biomass_costs_node[name].values * biomass_pot_node[name].values).sum() / biomass_pot_node[name].values.sum()).round(3)
        if snakemake.config['biomass']['nobiomass'] and name not in ['municipal solid waste']:
            biomass_potential[name] = 0


    print('municipal solid waste cost: ',biomass_costs['municipal solid waste'])
    n.madd("Generator",
           nodes + " municipal solid waste",
           bus=nodes + " municipal solid waste",
           carrier="municipal solid waste",
           p_nom=biomass_potential['municipal solid waste'] / 8760,
           marginal_cost=biomass_costs['municipal solid waste'])


    digestible_biomass_types = ["manureslurry", "sewage sludge", "straw"]
    for name in digestible_biomass_types:
        print(name,' cost: ',biomass_costs[name])

        n.add("Carrier", name + " digestible biomass")

        n.madd("Generator",
               nodes + " " + name + " digestible biomass",
               bus=nodes + " digestible biomass",
               carrier=name + " digestible biomass",
               p_nom_extendable=True,
               p_nom_max=biomass_potential[name] / 8760,
               marginal_cost=biomass_costs[name])

    # TODO: gas grid cost added for biogas processes in insert_gas_distribution_costs, but demands refining! Also add CO2 transport cost!
    n.madd("Link",
           nodes + " biogas",
           bus0=nodes + " digestible biomass",
           bus1="EU gas",
           bus3="co2 atmosphere",
           carrier="biogas",
           capital_cost=(costs.at["biogas", "fixed"] + costs.at["biogas upgrading", "fixed"]) * costs.at["biogas","efficiency"],
           marginal_cost=costs.at["biogas upgrading", "VOM"] * costs.at["biogas","efficiency"],
           efficiency=costs.at["biogas","efficiency"],
           efficiency3=-costs.at['gas', 'CO2 intensity'] * costs.at["biogas","efficiency"],
           p_nom_extendable=True)

    if beccs:
        n.madd("Link",
               nodes + " biogas CC",
               bus0=nodes + " digestible biomass",
               bus1="EU gas",
               bus2="co2 stored",
               bus3="co2 atmosphere",
               
               carrier="biogas CC",
               capital_cost=(costs.at["biogas CC", "fixed"] + costs.at["biogas upgrading", "fixed"]) * costs.at["biogas CC", "efficiency"]
                            + costs.at['biomass CHP capture', 'fixed'] * costs.at["biogas CC", "CO2 stored"],
               # Assuming that the CO2 from upgrading is pure, such as in amine scrubbing. I.e., with and without CC is
               # equivalent. Adding biomass CHP capture because biogas is often small-scale and decentral so further
               # from e.g. CO2 grid or buyers. This is a proxy for the added cost for e.g. a raw biogas pipeline to a central upgrading facility
               marginal_cost=(costs.at["biogas CC", "VOM"] + costs.at["biogas upgrading", "VOM"]) * costs.at["biogas","efficiency"],
               efficiency=costs.at["biogas CC", "efficiency"],
               efficiency2=costs.at["biogas CC", "CO2 stored"] * options["cc_fraction"],
               efficiency3=-costs.at['gas', 'CO2 intensity'] * costs.at["biogas CC", "efficiency"] - costs.at["biogas CC", "CO2 stored"] * options["cc_fraction"],
               
               p_nom_extendable=True)


    solid_biomass_types = ["forest residues", "industry wood residues", "landscape care"]
    for name in solid_biomass_types:
        n.add("Carrier", name + " solid biomass")

        # #Update solid biomass costs according to sensitivity settings
        # if 'BM0' in snakemake.wildcards.bm_s:
        #     pass
        # elif 'BM2' in snakemake.wildcards.bm_s:
        #     biomass_costs[name] = biomass_costs[name] + (biomass_import_price - biomass_costs[name]) / 3
        # elif 'BM3' in snakemake.wildcards.bm_s:
        #     biomass_costs[name] = biomass_costs[name] + (biomass_import_price - biomass_costs[name]) * 2 / 3
        # elif 'BM4' in snakemake.wildcards.bm_s:
        #     biomass_costs[name] = biomass_import_price
        # elif 'BM1' in snakemake.wildcards.bm_s:
        #     pass

        print(name, ' cost: ', biomass_costs[name])

        n.madd("Generator",
               nodes + " " + name + " solid biomass",
               bus=nodes + " solid biomass",
               carrier=name + " solid biomass",
               p_nom_extendable=True,
               p_nom_max=biomass_potential[name] / 8760,
               marginal_cost=biomass_costs[name])

    for o in opts:
        if o[o.find("B") + 4:o.find("B") + 6] == "Im" and not snakemake.config['biomass']['nobiomass']:
            if 'IM0' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0
            elif 'IM1' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0
            elif 'IM2' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0.01
            elif 'IM3' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0.05
            elif 'IM4' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0.1
            elif 'IM5' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0.15
            elif 'IM6' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0.2
            elif 'IM7' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0.25
            elif 'IM8' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0.35
            elif 'IM9' in snakemake.wildcards.bmim_s:
                bm_im_ef = 0.4
            elif 'IM10' in snakemake.wildcards.bmim_s:
                bm_im_ef = -0.05
            elif 'IM11' in snakemake.wildcards.bmim_s:
                bm_im_ef = -0.1
            elif 'IM12' in snakemake.wildcards.bmim_s:
                bm_im_ef = -0.15
            elif 'IM13' in snakemake.wildcards.bmim_s:
                bm_im_ef = -0.2

            print("Adding biomass import with cost ", biomass_import_price, ' EUR/MWh, and embedded emissions of ', bm_im_ef*100, '%')

            n.add("Carrier", "solid biomass import")

            n.madd("Bus",
                   ["solid biomass import"],
                   location="EU",
                   carrier="solid biomass import")

            import_potential = {}
            import_cost = {}
            import_name = {}
            superfluous = {}
            tot_EU_biomass = biomass_pot_node.values.sum() - biomass_pot_node["not included"].values.sum()
            print('Total EU biomass: ', tot_EU_biomass * 3.6 / 1e9, ' EJ')
            step_size = 1  # EJ
            biomass_import_limit_low_level = 20e9  # EJ

            for num in range(1, 20):
                import_name[num] = "import" + str(num)
                if num == 1:
                    import_potential[num] = max(biomass_import_limit_low_level / 3.6 - tot_EU_biomass,
                                                0)  # substract EU biomass from 20 EJ. If EU biomass > 20, return 0
                    import_cost[num] = biomass_import_price  # EUR/MWh
                    superfluous = min(biomass_import_limit_low_level / 3.6 - tot_EU_biomass,
                                      0)  # If EU biomass > 20, reduce the following group(s)
                else:
                    import_potential[num] = max(step_size * 1e9 / 3.6 + superfluous, 0)  # EJ --> MWh
                    import_cost[num] = biomass_import_price + (step_size * 0.25 * (num - 1)) * 3.6  # EUR/MWh
                    superfluous += min(-superfluous, step_size * 1e9 / 3.6)


                n.madd("Store",
                       [import_name[num] + " solid biomass"],
                       bus="solid biomass import",
                       e_nom_extendable=True,
                       e_cyclic=True,
                       carrier="solid biomass import",
                       )

                n.madd("Generator",
                       [import_name[num] + " solid biomass"],
                       bus="solid biomass import",
                       carrier="solid biomass import",
                       p_nom_extendable=True,
                       p_nom_max=import_potential[num] / 8760,
                       marginal_cost=import_cost[num])

                n.madd("Link",
                       nodes + " " + import_name[num] + " solid biomass",
                       bus0="solid biomass import",
                       bus1=nodes + " solid biomass",
                       bus3="co2 atmosphere",
                       carrier="solid biomass import",
                       efficiency=1.,
                       efficiency3=bm_im_ef * costs.at['solid biomass', 'CO2 intensity'],
                       p_nom_extendable=True)

    n.madd("Link",
           nodes + " solid biomass to gas",
           bus0=nodes + " solid biomass",
           bus1="EU gas",
           bus3="co2 atmosphere",
           carrier="BioSNG",
           lifetime=costs.at['BioSNG', 'lifetime'],
           efficiency=costs.at['BioSNG', 'efficiency'],
           efficiency3=-costs.at['solid biomass', 'CO2 intensity'] + costs.at['BioSNG', 'CO2 stored'],
           p_nom_extendable=True,
           capital_cost=costs.at['BioSNG', 'fixed'] * costs.at['BioSNG', 'efficiency'],
           marginal_cost=costs.at['BioSNG', 'efficiency']*costs.loc["BioSNG", "VOM"]
           )

    if beccs:
        n.madd("Link",
               nodes + " solid biomass to gas CC",
               bus0=nodes + " solid biomass",
               bus1="EU gas",
               bus2="co2 stored",
               bus3="co2 atmosphere",
               
               carrier="BioSNG CC",
               lifetime=costs.at['BioSNG', 'lifetime'],
               efficiency=costs.at['BioSNG', 'efficiency'],
               efficiency2=costs.at['BioSNG', 'CO2 stored'] * options["cc_fraction"],
               efficiency3=-costs.at['solid biomass', 'CO2 intensity'] + costs.at['BioSNG', 'CO2 stored'] * (1 - options["cc_fraction"]),
               
               p_nom_extendable=True,
               capital_cost=costs.at['BioSNG', 'fixed'] * costs.at['BioSNG', 'efficiency'] + costs.at['biomass CHP capture', 'fixed'] * costs.at[
                   "BioSNG", "CO2 stored"],
               marginal_cost=costs.at['BioSNG', 'efficiency']*costs.loc["BioSNG", "VOM"]
               )

        n.madd("Link",
               nodes + " solid biomass to hydrogen CC",
               bus0=nodes + " solid biomass",
               bus1=nodes + " H2",
               bus2="co2 stored",
               bus3="co2 atmosphere",
               
               carrier="solid biomass to hydrogen CC",
               efficiency=costs.at['solid biomass to hydrogen', 'efficiency'],
               efficiency2=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
               efficiency3=-costs.at['solid biomass', 'CO2 intensity'] + costs.at['solid biomass', 'CO2 intensity'] * (1 - options["cc_fraction"]),
               
               p_nom_extendable=True,
               capital_cost=costs.at['solid biomass to hydrogen', 'fixed'] * costs.at['solid biomass to hydrogen', 'efficiency']
                            + costs.at['biomass CHP capture', 'fixed'] * costs.at['solid biomass', 'CO2 intensity'],
               marginal_cost=0.,
               )

    n.madd("Link",
           nodes + " biomass to liquid",
           bus0=nodes + " solid biomass",
           bus1=spatial.oil.nodes,
           bus3="co2 atmosphere",
           carrier="biomass to liquid",
           lifetime=costs.at['BtL', 'lifetime'],
           efficiency=costs.at['BtL', 'efficiency'],
           efficiency3=-costs.at['solid biomass', 'CO2 intensity'] + costs.at['BtL', 'CO2 stored'],
           p_nom_extendable=True,
           capital_cost=costs.at['BtL', 'fixed'] * costs.at['BtL', 'efficiency'],
           marginal_cost=costs.at['BtL', 'efficiency']*costs.loc["BtL", "VOM"]
           )

    if beccs:
        #Assuming that acid gas removal (incl. CO2) from syngas i performed with Rectisol process (Methanol) and electricity demand for this is included in the base process
        n.madd("Link",
               nodes + " biomass to liquid CC",
               bus0=nodes + " solid biomass",
               bus1=spatial.oil.nodes,
               bus2="co2 stored",
               bus3="co2 atmosphere",
               carrier="biomass to liquid CC",
               lifetime=costs.at['BtL', 'lifetime'],
               efficiency=costs.at['BtL', 'efficiency'],
               efficiency2=costs.at['BtL', 'CO2 stored'] * options["cc_fraction"],
               efficiency3=-costs.at['solid biomass', 'CO2 intensity'] + costs.at['BtL', 'CO2 stored'] * (1 - options["cc_fraction"]),
               p_nom_extendable=True,
               capital_cost=costs.at['BtL', 'fixed'] * costs.at['BtL', 'efficiency'] + costs.at['biomass CHP capture', 'fixed'] * costs.at[
                   "BtL", "CO2 stored"],
               marginal_cost=costs.at['BtL', 'efficiency'] * costs.loc["BtL", "VOM"]
               )

    if options['electrobiofuels']:
        print('Adding electrobiofuels')
        efuel_scale_factor = costs.at['BtL', 'C stored']
        n.madd("Link",
               nodes + " electrobiofuels",
               bus0=nodes + " solid biomass",
               bus1=spatial.oil.nodes,
               bus5=nodes + " H2",
               bus3="co2 atmosphere",
               carrier="electrobiofuels",
               lifetime=costs.at['electrobiofuels', 'lifetime'],
               efficiency=costs.at['electrobiofuels', 'efficiency-biomass'],
               efficiency5=-costs.at['electrobiofuels', 'efficiency-hydrogen'],
               efficiency3=-costs.at['solid biomass', 'CO2 intensity'] + costs.at['BtL', 'CO2 stored'] * (1 - costs.at['Fischer-Tropsch', 'capture rate']),
               p_nom_extendable=True,
               capital_cost=costs.at['BtL', 'fixed'] * costs.at['BtL', 'efficiency'] / costs.at['electrobiofuels', 'efficiency-biomass'] \
                            + efuel_scale_factor * costs.at['Fischer-Tropsch', 'fixed'] * costs.at['Fischer-Tropsch', 'efficiency'] / costs.at['electrobiofuels', 'efficiency-hydrogen'],
               marginal_cost=costs.at['BtL', 'VOM'] * costs.at['BtL', 'efficiency'] / costs.at['electrobiofuels', 'efficiency-biomass'] \
                             + efuel_scale_factor * costs.at['Fischer-Tropsch', 'VOM'] * costs.at['Fischer-Tropsch', 'efficiency'] / costs.at['electrobiofuels', 'efficiency-hydrogen']
               )


    # TODO: Add real data for bioelectricity without CHP!
    n.madd("Link",
           nodes + " solid biomass to electricity",
           bus0=nodes + " solid biomass",
           bus1=nodes,
           bus3="co2 atmosphere",
           carrier="solid biomass to electricity",
           p_nom_extendable=True,
           capital_cost=0.7 * costs.at['central solid biomass CHP', 'fixed'] * costs.at[
               'central solid biomass CHP', 'efficiency'],
           marginal_cost=costs.at['central solid biomass CHP', 'VOM'],
           efficiency=1.3 * costs.at['central solid biomass CHP', 'efficiency'],
           efficiency3=costs.at['solid biomass', 'CO2 intensity']-costs.at['solid biomass', 'CO2 intensity'],
           lifetime=costs.at['central solid biomass CHP', 'lifetime'])

    if beccs:
        
            
        n.madd("Link",
               nodes + " solid biomass to electricity CC",
               bus0=nodes + " solid biomass",
               bus1=nodes,
               bus2="co2 stored",
               bus3="co2 atmosphere",
               bus4= "EU MEA",
               carrier="solid biomass to electricity CC",
               p_nom_extendable=True,
               capital_cost=0.7 * costs.at['central solid biomass CHP CC', 'fixed'] * costs.at[
                   'central solid biomass CHP CC', 'efficiency']
                            + costs.at['biomass CHP capture', 'fixed'] * costs.at['solid biomass', 'CO2 intensity'],
               marginal_cost=costs.at['central solid biomass CHP CC', 'VOM'],
               efficiency=1.3 * costs.at['central solid biomass CHP CC', 'efficiency'],
               efficiency2=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
               efficiency3=-costs.at['solid biomass', 'CO2 intensity'] + costs.at['solid biomass', 'CO2 intensity'] * (1 - options["cc_fraction"]),
               efficiency4= - costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"], #3.98e-3,
               # p_nom_ratio=costs.at['central solid biomass CHP', 'p_nom_ratio'],
               lifetime=costs.at['central solid biomass CHP CC', 'lifetime'])

    # AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty and options["chp"]:
        urban_central = urban_central.str[:-len(" urban central heat")]

        n.madd("Link",
               urban_central + " urban central solid biomass CHP",
               bus0=urban_central + " solid biomass",
               bus1=urban_central,
               bus4=urban_central + " urban central heat",
               bus3="co2 atmosphere",
               carrier="urban central solid biomass CHP",
               p_nom_extendable=True,
               capital_cost=costs.at['central solid biomass CHP', 'fixed'] * costs.at[
                   'central solid biomass CHP', 'efficiency'],
               marginal_cost=costs.at['central solid biomass CHP', 'VOM'],
               efficiency=costs.at['central solid biomass CHP', 'efficiency'],
               efficiency4=costs.at['central solid biomass CHP', 'efficiency-heat'],
               efficiency3=costs.at['solid biomass', 'CO2 intensity']-costs.at['solid biomass', 'CO2 intensity'],
               c_b=costs.at['central solid biomass CHP', 'c_b'],
               c_v=costs.at['central solid biomass CHP', 'c_v'],
               lifetime=costs.at['central solid biomass CHP', 'lifetime'])

        n.madd("Link",
               urban_central + " urban central solid biomass heat",
               bus0=urban_central + " solid biomass",
               bus1=urban_central + " urban central heat",
               carrier="solid biomass district heat",
               p_nom_extendable=True,
               efficiency=costs.at['solid biomass boiler steam', 'efficiency'],
               capital_cost=costs.at['solid biomass boiler steam', 'fixed'] * costs.at[
                   'solid biomass boiler steam', 'efficiency'],
               marginal_cost=costs.at['solid biomass boiler steam', 'VOM'])

        if beccs:
        
            
                
            n.madd("Link",
                   urban_central + " urban central solid biomass CHP CC",
                   bus0=urban_central + " solid biomass",
                   bus1=urban_central,
                   bus2="co2 stored",
                   bus3="co2 atmosphere",
                   bus4=urban_central + " urban central heat",
                   bus5="EU MEA",
                   carrier="urban central solid biomass CHP CC",
                   p_nom_extendable=True,
                   capital_cost=costs.at['central solid biomass CHP CC', 'fixed'] * costs.at['central solid biomass CHP CC', 'efficiency']
                                + costs.at['biomass CHP capture', 'fixed'] * costs.at['solid biomass', 'CO2 intensity'],
                   marginal_cost=costs.at['central solid biomass CHP CC', 'VOM'],
                   efficiency=costs.at['central solid biomass CHP CC', 'efficiency'],
                   efficiency2=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
                   efficiency3=costs.at['solid biomass', 'CO2 intensity'] * (1 - options["cc_fraction"])-costs.at['solid biomass', 'CO2 intensity'],
                   efficiency4=costs.at['central solid biomass CHP CC', 'efficiency-heat'],
                   efficiency5= - costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"]* options["MEA_consumption"], #3.98e-3,
                   c_b=costs.at['central solid biomass CHP CC', 'c_b'],
                   c_v=costs.at['central solid biomass CHP CC', 'c_v'],
                   lifetime=costs.at['central solid biomass CHP CC', 'lifetime'])

            if options['solid_biomass_powerboost_CC']:
                n.madd("Link",
                       urban_central + " urban central solid biomass CHP powerboost CC",
                       bus0=urban_central + " solid biomass",
                       bus1=urban_central,
                       bus4=urban_central + " urban central heat",
                       bus3="co2 atmosphere",
                       bus2="co2 stored",
                       carrier="urban central solid biomass CHP powerboost CC",
                       p_nom_extendable=True,
                       capital_cost=costs.at['central solid biomass CHP powerboost CC', 'fixed'] * costs.at['central solid biomass CHP powerboost CC', 'efficiency']
                                    + costs.at['biomass CHP capture', 'fixed'] * costs.at['solid biomass', 'CO2 intensity'],
                       marginal_cost=costs.at['central solid biomass CHP powerboost CC', 'VOM'],
                       efficiency=costs.at['central solid biomass CHP powerboost CC', 'efficiency'],
                       efficiency4=costs.at['central solid biomass CHP powerboost CC', 'efficiency-heat'],
                       efficiency3=costs.at['solid biomass', 'CO2 intensity'] * (1 - options["cc_fraction"])-costs.at['solid biomass', 'CO2 intensity'],
                       efficiency2=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
                       c_b=costs.at['central solid biomass CHP powerboost CC', 'c_b'],
                       c_v=costs.at['central solid biomass CHP powerboost CC', 'c_v'],
                       lifetime=costs.at['central solid biomass CHP powerboost CC', 'lifetime'])

            n.madd("Link",
                   urban_central + " urban central solid biomass heat CC",
                   bus0=urban_central + " solid biomass",
                   bus1=urban_central + " urban central heat",
                   bus3="co2 atmosphere",
                   bus2="co2 stored",
                   bus4="EU MEA",
                   carrier="solid biomass district heat CC",
                   p_nom_extendable=True,
                   efficiency=costs.at['solid biomass boiler steam CC', 'efficiency'],
                   capital_cost=costs.at['solid biomass boiler steam CC', 'fixed'] * costs.at[
                       'solid biomass boiler steam CC', 'efficiency'] + costs.at[
                                    "biomass CHP capture", "fixed"] * costs.at['solid biomass', 'CO2 intensity'],
                   marginal_cost=costs.at['solid biomass boiler steam CC', 'VOM'],
                   efficiency3=costs.at['solid biomass', 'CO2 intensity'] * (
                           1 - options["cc_fraction"]) - costs.at[
                                   'solid biomass', 'CO2 intensity'],
                   efficiency2=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
                   efficiency4= -costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"], #3.98e-3,
                   lifetime=costs.at['solid biomass boiler steam CC', 'lifetime'])

        if options['waste_chp']:
            print('Adding waste CHPs')
            n.madd("Link",
                   urban_central + " waste CHP",
                   bus0=urban_central + " municipal solid waste",
                   bus1=urban_central,
                   bus4=urban_central + " urban central heat",
                   bus3="co2 atmosphere",
                   carrier="urban central waste incineration",
                   p_nom_extendable=True,
                   # p_nom=biomass_potential['municipal solid waste'] / 8760,
                   capital_cost=costs.at['waste CHP', 'fixed'] * costs.at['waste CHP', 'efficiency'],
                   marginal_cost=costs.at['waste CHP', 'VOM'],
                   efficiency=costs.at['waste CHP', 'efficiency'],
                   efficiency4=costs.at['waste CHP', 'efficiency-heat'],
                   efficiency3=costs.at['solid biomass', 'CO2 intensity']-costs.at['solid biomass', 'CO2 intensity'],
                   lifetime=costs.at['waste CHP', 'lifetime'])

            if beccs:
                n.madd("Link",
                       urban_central + " waste CHP CC",
                       bus0=urban_central + " municipal solid waste",
                       bus1=urban_central,
                       bus4=urban_central + " urban central heat",
                       bus3="co2 atmosphere",
                       bus2="co2 stored",
                       bus5="EU MEA",
                       carrier="urban central waste incineration CC",
                       p_nom_extendable=True,
                       # p_nom=costs.at['waste CHP CC', 'efficiency'] * biomass_potential['municipal solid waste'] / 8760,
                       capital_cost=costs.at['waste CHP CC', 'fixed'] * costs.at['waste CHP CC', 'efficiency']
                                    + costs.at['biomass CHP capture', 'fixed'] * costs.at['solid biomass', 'CO2 intensity'],
                       marginal_cost=costs.at['waste CHP CC', 'VOM'],
                       efficiency=costs.at['waste CHP CC', 'efficiency'],
                       efficiency4=costs.at['waste CHP CC', 'efficiency-heat'],
                       #Assuming same CO2 intensity as solid biomass
                       efficiency3=costs.at['solid biomass', 'CO2 intensity'] * (1 - options["cc_fraction"])-costs.at['solid biomass', 'CO2 intensity'],
                       efficiency2=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
                       efficiency5= - costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"], #3.98e-3
                       c_b=costs.at['waste CHP CC', 'c_b'],
                       c_v=costs.at['waste CHP CC', 'c_v'],
                       lifetime=costs.at['waste CHP CC', 'lifetime'])

    if options["biomass_boiler"]:
        nodes_heat = create_nodes_for_heat_sector()[0]
        for name in ["residential rural", "services rural",
                     "residential urban decentral", "services urban decentral"]:

            n.madd("Link",
                nodes_heat[name] + f" {name} biomass boiler",
                p_nom_extendable=True,
                bus0=spatial.biomass.df.loc[nodes_heat[name], "nodes"].values,
                bus1=nodes_heat[name] + f" {name} heat",
                carrier=name + " biomass boiler",
                efficiency=costs.at['biomass boiler', 'efficiency'],
                capital_cost=costs.at['biomass boiler', 'efficiency'] * costs.at['biomass boiler', 'fixed'],
                marginal_cost=costs.at['biomass boiler', 'pelletizing cost'],
                lifetime=costs.at['biomass boiler', 'lifetime']
            )


def add_industry(n, costs):

    logger.info("Add industrial demand")
    nodes = pop_layout.index

    # 1e6 to convert TWh to MWh
    industrial_demand = pd.read_csv(snakemake.input.industrial_demand, index_col=0) * 1e6

    n.madd("Bus",
           nodes + " >120C process steam",
           location=nodes,
           carrier=">120C process steam",
           unit="MWh_LHV")

    # if options["biomass_transport"]:
    #     p_set = industrial_demand.loc[spatial.biomass.locations, "solid biomass"].rename(index=lambda x: x + " solid biomass for industry") / 8760
    # else:
    #     p_set = industrial_demand["solid biomass"].sum() / 8760

    n.madd("Load",
           nodes,
           suffix=" >120C process steam",
           bus=nodes + " >120C process steam",
           carrier=">120C process steam",
           p_set=industrial_demand.loc[nodes, "solid biomass"] / 8760.)

    n.madd("Bus",
           nodes + " mediumT industry",
           location=nodes,
           carrier="mediumT industry",
           unit="MWh_LHV")

    #TODO: Set real shares of medium and high T industry
    n.madd("Load",
           nodes,
           suffix=" mediumT industry",
           bus=nodes + " mediumT industry",
           carrier="mediumT industry",
           p_set=0.36*industrial_demand.loc[nodes, "methane"] / 8760.)

    n.madd("Bus",
           nodes + " highT industry",
           location=nodes,
           carrier="highT industry")

    n.madd("Load",
           nodes,
           suffix=" highT industry",
           bus=nodes + " highT industry",
           carrier="highT industry",
           p_set=0.64*industrial_demand.loc[nodes, "methane"] / 8760.)

    for o in opts:
        if "B" in o:
            
            if snakemake.config['biomass']['>120C industry steam biomass']:
                n.madd("Link",
                       nodes,
                       suffix=" solid biomass for >120C industry",
                       bus0=nodes + " solid biomass",
                       bus1=nodes + " >120C process steam",
                       bus3="co2 atmosphere",
                       carrier=">120C process steam solid biomass",
                       p_nom_extendable=True,
                       p_min_pu=0.8,
                       efficiency=costs.at['solid biomass boiler steam', 'efficiency'],
                       efficiency3=costs.at['solid biomass', 'CO2 intensity']-costs.at['solid biomass', 'CO2 intensity'],
                       capital_cost=costs.at['solid biomass boiler steam', 'fixed'] * costs.at['solid biomass boiler steam', 'efficiency'],
                       marginal_cost=costs.at['solid biomass boiler steam', 'VOM'])

            if snakemake.config['biomass']['mediumT industry biomass']:
                n.madd("Link",
                       nodes,
                       suffix=" solid biomass for mediumT industry",
                       bus0=nodes + " solid biomass",
                       bus1=nodes + " mediumT industry",
                       bus3="co2 atmosphere",
                       carrier="solid biomass for mediumT industry",
                       p_nom_extendable=True,
                       p_min_pu=0.8,
                       efficiency=costs.at['direct firing solid fuels', 'efficiency'],
                       efficiency3=costs.at['solid biomass', 'CO2 intensity']-costs.at['solid biomass', 'CO2 intensity'],
                       capital_cost=costs.at['direct firing solid fuels', 'fixed'] * costs.at['direct firing solid fuels', 'efficiency'],
                       marginal_cost=costs.at['direct firing solid fuels', 'VOM'] + costs.at['biomass boiler', 'pelletizing cost'],
                       lifetime=costs.at['direct firing solid fuels', 'lifetime'])


            if snakemake.config['biomass']['beccs']:
                if snakemake.config['biomass']['>120C industry steam biomass']:
                    n.madd("Link",
                           nodes,
                           suffix=" solid biomass for >120C industry CC",
                           bus0=nodes + " solid biomass",
                           bus1=nodes + " >120C process steam",
                           bus3="co2 atmosphere",
                           bus2="co2 stored",
                           bus4="EU MEA",
                           carrier=">120C process steam solid biomass CC",
                           p_nom_extendable=True,
                           p_min_pu=0.8,
                           efficiency=costs.at['solid biomass boiler steam CC', 'efficiency'],
                           capital_cost=costs.at['solid biomass boiler steam CC', 'fixed'] * costs.at['solid biomass boiler steam CC', 'efficiency'] + costs.at[
                               "biomass CHP capture", "fixed"] * costs.at['solid biomass', 'CO2 intensity'],
                           marginal_cost=costs.at['solid biomass boiler steam CC', 'VOM'],
                           efficiency3=costs.at['solid biomass', 'CO2 intensity'] * (
                                   1 - options["cc_fraction"])-costs.at['solid biomass', 'CO2 intensity'],
                           efficiency2=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
                           efficiency4= - costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"],#3.98e-3,
                           lifetime=costs.at['solid biomass boiler steam CC', 'lifetime'])

                if snakemake.config['biomass']['mediumT industry biomass']:
                    n.madd("Link",
                           nodes,
                           suffix=" solid biomass for mediumT industry CC",
                           bus0=nodes + " solid biomass",
                           bus1=nodes + " mediumT industry",
                           bus3="co2 atmosphere",
                           bus2="co2 stored",
                           bus4="EU MEA",
                           carrier="solid biomass for mediumT industry CC",
                           p_nom_extendable=True,
                           p_min_pu=0.8,
                           efficiency=costs.at['direct firing solid fuels CC', 'efficiency'],
                           capital_cost=costs.at['direct firing solid fuels CC', 'fixed'] * costs.at['direct firing solid fuels CC', 'efficiency'] + costs.at[
                               "biomass CHP capture", "fixed"] * costs.at['solid biomass', 'CO2 intensity'],
                           marginal_cost=costs.at['direct firing solid fuels CC', 'VOM'] + costs.at['biomass boiler', 'pelletizing cost'],
                           efficiency3=costs.at['solid biomass', 'CO2 intensity'] * (
                                   1 - options["cc_fraction"])-costs.at['solid biomass', 'CO2 intensity'],
                           efficiency2=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
                           efficiency4= - costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"], #3.98e-3,
                           lifetime=costs.at['direct firing solid fuels CC', 'lifetime'])

    if options["industrial_steam_methane"]:
    
        
             
        n.madd("Link",
               nodes,
               suffix=" gas for >120C industry",
               bus0="EU gas",
               bus1=nodes + " >120C process steam",
               bus3="co2 atmosphere",
               carrier=">120C process steam methane",
               p_nom_extendable=True,
               p_min_pu=0.8,
               efficiency=costs.at['gas boiler steam', 'efficiency'],
               capital_cost=costs.at['gas boiler steam', 'fixed'] * costs.at['gas boiler steam', 'efficiency'],
               marginal_cost=costs.at['gas boiler steam', 'VOM'],
               efficiency3=costs.at['gas', 'CO2 intensity'])

        eta = costs.at['gas boiler steam', 'efficiency'] - costs.at['gas', 'CO2 intensity'] * costs.at['biomass CHP capture', 'heat-input']
        n.madd("Link",
               nodes,
               suffix=" gas for >120C industry CC",
               bus0="EU gas",
               bus1=nodes + " >120C process steam",
               bus3="co2 atmosphere",
               bus2="co2 stored",
               bus4="EU MEA",
               carrier=">120C process steam methane CC",
               p_nom_extendable=True,
               p_min_pu=0.8,
               efficiency=eta,
               capital_cost=costs.at['gas boiler steam', 'fixed'] * costs.at['gas boiler steam', 'efficiency'] + costs.at["biomass CHP capture", "fixed"] * costs.at[
                   'gas', 'CO2 intensity'],
               marginal_cost=costs.at['gas boiler steam', 'VOM'],
               efficiency3=costs.at['gas', 'CO2 intensity'] * (1 - options["cc_fraction"]),
               efficiency2=costs.at['gas', 'CO2 intensity'] * options["cc_fraction"],
               efficiency4=-costs.at['gas', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"], #costs.at["MEA consumption", "efficiency"], #3.98e-3,
               lifetime=costs.at['gas boiler steam', 'lifetime'])

    if options["industrial_steam_heat_pumps"]:

        n.madd("Bus",
               nodes + " <120C process steam",
               location=nodes,
               carrier="<120C process steam",
               unit="MWh_LHV")

        eta=costs.at['industrial heat pump medium temperature', 'efficiency'] * options['industrial_steam_heat_pump_factor']
        print('adding industrial heat pump with COP: ', eta)

        n.madd("Link",
               nodes,
               suffix=" industrial heat pump steam for <120C industry",
               bus0=nodes,
               bus1=nodes + " <120C process steam",
               carrier="<120C process steam heat pump",
               p_nom_extendable=True,
               p_min_pu=0.8,
               efficiency=eta,
               capital_cost=costs.at['industrial heat pump medium temperature', 'fixed'] * eta,
               marginal_cost=costs.at['industrial heat pump medium temperature', 'VOM'],
               lifetime=costs.at['industrial heat pump medium temperature', 'lifetime'])

    if options["industrial_steam_electric_boiler"]:
        n.madd("Link",
               nodes,
               suffix=" electricity for >120C industry",
               bus0=nodes,
               bus1=nodes + " >120C process steam",
               carrier=">120C process steam electricity",
               p_nom_extendable=True,
               p_min_pu=0.8,
               efficiency=costs.at['electric boiler steam', 'efficiency'],
               capital_cost=costs.at['electric boiler steam', 'fixed'] * costs.at['electric boiler steam', 'efficiency'],
               marginal_cost=costs.at['electric boiler steam', 'VOM'],
               lifetime=costs.at['electric boiler steam', 'lifetime'])

    if options["gas_for_mediumT_industry"]:
        #TODO: add electricity input from DEA and adapt VOM to exclude electricity cost!
        n.madd("Link",
               nodes,
               suffix=" gas for mediumT industry",
               bus0="EU gas",
               bus1=nodes + " mediumT industry",
               bus3="co2 atmosphere",
               carrier="gas for mediumT industry",
               p_nom_extendable=True,
               p_min_pu=0.8,
               efficiency=costs.at['direct firing gas', 'efficiency'],
               efficiency3=costs.at['gas', 'CO2 intensity'],
               capital_cost=costs.at['direct firing gas', 'fixed'] * costs.at['direct firing gas', 'efficiency'],
               marginal_cost=costs.at['direct firing gas', 'VOM'],
               lifetime=costs.at['direct firing gas', 'lifetime'])

        eta=costs.at['direct firing gas', 'efficiency'] - costs.at['gas', 'CO2 intensity'] * costs.at['biomass CHP capture', 'heat-input']
        
        
            
        n.madd("Link",
               nodes,
               suffix=" gas for mediumT industry CC",
               bus0="EU gas",
               bus1=nodes + " mediumT industry",
               bus3="co2 atmosphere",
               bus2="co2 stored",
               bus4="EU MEA",
               carrier="gas for mediumT industry CC",
               p_nom_extendable=True,
               p_min_pu=0.8,
               efficiency=eta,
               efficiency3=costs.at['gas', 'CO2 intensity'] * (1 - options["cc_fraction"]),
               efficiency2=costs.at['gas', 'CO2 intensity'] * options["cc_fraction"],
               efficiency4= - costs.at['gas', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"],#3.98e-3,
               capital_cost=costs.at['direct firing gas CC', 'fixed'] * costs.at['direct firing gas CC', 'efficiency'] + costs.at["biomass CHP capture", "fixed"] * costs.at['gas', 'CO2 intensity'],
               marginal_cost=costs.at['direct firing gas CC', 'VOM'],
               lifetime=costs.at['direct firing gas', 'lifetime'])

    #TODO: differentiate medium and high T processes
    n.madd("Link",
           nodes,
           suffix=" gas for highT industry",
           bus0="EU gas",
           bus1=nodes + " highT industry",
           bus3="co2 atmosphere",
           carrier="gas for highT industry",
           p_nom_extendable=True,
           p_min_pu=0.8,
           efficiency=costs.at['direct firing gas', 'efficiency'],
           efficiency3=costs.at['gas', 'CO2 intensity'],
           capital_cost=costs.at['direct firing gas', 'fixed'] * costs.at['direct firing gas', 'efficiency'],
           marginal_cost=costs.at['direct firing gas', 'VOM'],
           lifetime=costs.at['direct firing gas', 'lifetime'])

    eta=costs.at['direct firing gas', 'efficiency'] - costs.at['gas', 'CO2 intensity'] * costs.at['biomass CHP capture', 'heat-input']
    
    n.madd("Link",
           nodes,
           suffix=" gas for highT industry CC",
           bus0="EU gas",
           bus1=nodes + " highT industry",
           bus3="co2 atmosphere",
           bus2="co2 stored",
           bus4 ="EU MEA",
           carrier="gas for highT industry CC",
           p_nom_extendable=True,
           p_min_pu=0.8,
           efficiency=eta,
           efficiency3=costs.at['gas', 'CO2 intensity'] * (1 - options["cc_fraction"]),
           efficiency2=costs.at['gas', 'CO2 intensity'] * options["cc_fraction"],
           efficiency4= - costs.at['gas', 'CO2 intensity'] * options["cc_fraction"] * options["MEA_consumption"], #3.98e-3,
           capital_cost=costs.at['direct firing gas CC', 'fixed'] * costs.at['direct firing gas CC', 'efficiency'] + costs.at["biomass CHP capture", "fixed"] * costs.at['gas', 'CO2 intensity'],
           
           marginal_cost=costs.at['direct firing gas CC', 'VOM'],
           lifetime=costs.at['direct firing gas', 'lifetime'])

    if options["hydrogen_for_mediumT_industry"]:
        print('Adding H2 for mediumT industry')
        n.madd("Link",
               nodes,
               suffix=" hydrogen for mediumT industry",
               bus0=nodes + " H2",
               bus1=nodes + " mediumT industry",
               carrier="hydrogen for mediumT industry",
               capital_cost=10 * costs.at['direct firing gas', 'fixed'] * costs.at['direct firing gas', 'efficiency'],
               marginal_cost=10 * costs.at['direct firing gas', 'VOM'],
               p_nom_extendable=True,
               p_min_pu=0.8,
               efficiency=costs.at['direct firing gas', 'efficiency'])

    if options["hydrogen_for_highT_industry"]:
        print('Adding H2 for highT industry')
        n.madd("Link",
               nodes,
               suffix=" hydrogen for highT industry",
               bus0=nodes + " H2",
               bus1=nodes + " highT industry",
               carrier="hydrogen for highT industry",
               capital_cost=10 * costs.at['direct firing gas', 'fixed'] * costs.at['direct firing gas', 'efficiency'],
               marginal_cost=10 * costs.at['direct firing gas', 'VOM'],
               p_nom_extendable=True,
               p_min_pu=0.8,
               efficiency=costs.at['direct firing gas', 'efficiency'])

    n.madd("Load",
        nodes,
        suffix=" H2 for industry",
        bus=nodes + " H2",
        carrier="H2 for industry",
        p_set=industrial_demand.loc[nodes, "hydrogen"] / 8760
    )

    if options["shipping_hydrogen_liquefaction"]:

        n.madd("Bus",
            nodes,
            suffix=" H2 liquid",
            carrier="H2 liquid",
            location=nodes,
            unit="MWh_LHV"
        )

        n.madd("Link",
            nodes + " H2 liquefaction",
            bus0=nodes + " H2",
            bus1=nodes + " H2 liquid",
            carrier="H2 liquefaction",
            efficiency=costs.at["H2 liquefaction", 'efficiency'],
            capital_cost=costs.at["H2 liquefaction", 'fixed'],
            p_nom_extendable=True,
            lifetime=costs.at['H2 liquefaction', 'lifetime']
        )

        shipping_bus = nodes + " H2 liquid"
    else:
        shipping_bus = nodes + " H2"

    all_navigation = ["total international navigation", "total domestic navigation"]
    efficiency = options['shipping_average_efficiency'] / costs.at["fuel cell", "efficiency"]
    shipping_hydrogen_share = get(options['shipping_hydrogen_share'], investment_year)
    p_set = shipping_hydrogen_share * pop_weighted_energy_totals.loc[nodes, all_navigation].sum(axis=1) * 1e6 * efficiency / 8760

    n.madd("Load",
        nodes,
        suffix=" H2 for shipping",
        bus=shipping_bus,
        carrier="H2 for shipping",
        p_set=p_set
    )

    if shipping_hydrogen_share < 1:

        shipping_oil_share = 1 - shipping_hydrogen_share

        p_set = shipping_oil_share * get(options["shipping_demand"], investment_year) * pop_weighted_energy_totals.loc[nodes, all_navigation].sum(axis=1) * 1e6 / 8760.

        n.madd("Load",
            nodes,
            suffix=" shipping oil",
            bus=spatial.oil.nodes,
            carrier="shipping oil",
            p_set=p_set
        )

        co2 = shipping_oil_share * get(options["shipping_demand"], investment_year) * pop_weighted_energy_totals.loc[nodes, all_navigation].sum().sum() * 1e6 / 8760 * costs.at["oil", "CO2 intensity"]

        n.add("Load",
            "shipping oil emissions",
            bus="co2 atmosphere",
            carrier="shipping oil emissions",
            p_set=-co2
        )

    if "oil" not in n.buses.carrier.unique():
        n.madd("Bus",
            spatial.oil.nodes,
            location=spatial.oil.locations,
            carrier="oil",
            unit="MWh_LHV"
        )

    if "oil" not in n.stores.carrier.unique():

        #could correct to e.g. 0.001 EUR/kWh * annuity and O&M
        n.madd("Store",
            [oil_bus + " Store" for oil_bus in spatial.oil.nodes],
            bus=spatial.oil.nodes,
            e_nom_extendable=True,
            e_cyclic=True,
            carrier="oil",
        )

    if 'nofossil' in opts:
        extendable=False
    elif "oil" not in n.generators.carrier.unique():
        extendable=True

    if "oil" not in n.generators.carrier.unique():
        n.madd("Generator",
            spatial.oil.nodes,
            bus=spatial.oil.nodes,
            p_nom_extendable=extendable,
            carrier="oil",
            marginal_cost=costs.at["oil", 'fuel']
        )

    if options["oil_boilers"]:

        nodes_heat = create_nodes_for_heat_sector()[0]

        for name in ["residential rural", "services rural", "residential urban decentral", "services urban decentral"]:

            n.madd("Link",
                nodes_heat[name] + f" {name} oil boiler",
                p_nom_extendable=True,
                bus0=spatial.oil.nodes,
                bus1=nodes_heat[name] + f" {name}  heat",
                bus3="co2 atmosphere",
                carrier=f"{name} oil boiler",
                efficiency=costs.at['decentral oil boiler', 'efficiency'],
                efficiency3=costs.at['oil', 'CO2 intensity'],
                capital_cost=costs.at['decentral oil boiler', 'efficiency'] * costs.at['decentral oil boiler', 'fixed'],
                lifetime=costs.at['decentral oil boiler', 'lifetime']
            )

    n.madd("Link",
        nodes + " Fischer-Tropsch",
        bus0=nodes + " H2",
        bus1=spatial.oil.nodes,
        bus2=spatial.co2.nodes,
        bus3="co2 atmosphere",
        carrier="electrofuel",
        efficiency=costs.at["Fischer-Tropsch", 'efficiency'],
        capital_cost=costs.at["Fischer-Tropsch", 'fixed'] * costs.at["Fischer-Tropsch", 'efficiency'],
        efficiency2=-(2 - costs.at["Fischer-Tropsch", 'capture rate']) * costs.at['oil', 'CO2 intensity'] *
                       costs.at["Fischer-Tropsch", 'efficiency'],
        efficiency3=(1 - costs.at["Fischer-Tropsch", 'capture rate']) * costs.at['oil', 'CO2 intensity'] * costs.at[
               "Fischer-Tropsch", 'efficiency'],
        p_nom_extendable=True,
        lifetime=costs.at['Fischer-Tropsch', 'lifetime']
    )

    n.madd("Load",
        ["naphtha for industry"],
        bus=spatial.oil.nodes,
        carrier="naphtha for industry",
        p_set=industrial_demand.loc[nodes, "naphtha"].sum() / 8760
    )

    all_aviation = ["total international aviation", "total domestic aviation"]
    p_set = get(options["aviation_demand"], investment_year) * pop_weighted_energy_totals.loc[nodes, all_aviation].sum(axis=1).sum() * 1e6 / 8760

    n.madd("Load",
        ["kerosene for aviation"],
        bus=spatial.oil.nodes,
        carrier="kerosene for aviation",
        p_set=p_set
    )

    #NB: CO2 gets released again to atmosphere when plastics decay or kerosene is burned
    #except for the process emissions when naphtha is used for petrochemicals, which can be captured with other industry process emissions
    #tco2 per hour
    co2_release = ["naphtha for industry", "kerosene for aviation"]
    co2 = n.loads.loc[co2_release, "p_set"].sum() * costs.at["oil", 'CO2 intensity'] - industrial_demand.loc[nodes, "process emission from feedstock"].sum() / 8760

    n.add("Load",
        "oil emissions",
        bus="co2 atmosphere",
        carrier="oil emissions",
        p_set=-co2
    )

    # TODO simplify bus expression
    n.madd("Load",
        nodes,
        suffix=" low-temperature heat for industry",
        bus=[node + " urban central heat" if node + " urban central heat" in n.buses.index else node + " services urban decentral heat" for node in nodes],
        carrier="low-temperature heat for industry",
        p_set=industrial_demand.loc[nodes, "low-temperature heat"] / 8760
    )

    # remove today's industrial electricity demand by scaling down total electricity demand
    for ct in n.buses.country.dropna().unique():
        # TODO map onto n.bus.country
        loads_i = n.loads.index[(n.loads.index.str[:2] == ct) & (n.loads.carrier == "electricity")]
        if n.loads_t.p_set[loads_i].empty: continue
        factor = 1 - industrial_demand.loc[loads_i, "current electricity"].sum() / n.loads_t.p_set[loads_i].sum().sum()
        n.loads_t.p_set[loads_i] *= factor

    n.madd("Load",
        nodes,
        suffix=" industry electricity",
        bus=nodes,
        carrier="industry electricity",
        p_set=industrial_demand.loc[nodes, "electricity"] / 8760
    )

    n.add("Bus",
        "process emissions",
        location="EU",
        carrier="process emissions",
        unit="t_co2"
    )

    # this should be process emissions fossil+feedstock
    # then need load on atmosphere for feedstock emissions that are currently going to atmosphere via Link Fischer-Tropsch demand
    n.add("Load",
        "process emissions",
        bus="process emissions",
        carrier="process emissions",
        p_set=-industrial_demand.loc[nodes,["process emission", "process emission from feedstock"]].sum(axis=1).sum() / 8760
    )

    n.add("Link",
        "process emissions",
        bus0="process emissions",
        bus1="co2 atmosphere",
        carrier="process emissions",
        p_nom_extendable=True,
        efficiency=1.
    )

    #assume enough local waste heat for CC
    n.madd("Link",
        spatial.co2.locations,
        suffix=" process emissions CC",
        bus0="process emissions",
        bus1="co2 atmosphere",
        bus2=spatial.co2.nodes,
        bus3="EU MEA",
        carrier="process emissions CC",
        p_nom_extendable=True,
        capital_cost=costs.at["cement capture", "fixed"],
        efficiency=1 - options["cc_fraction"],
        efficiency2=options["cc_fraction"],
        efficiency3=-options["MEA_consumption"],
        lifetime=costs.at['cement capture', 'lifetime']
    )

    if options.get("ammonia"):

        if options["ammonia"] == 'regional':
            p_set = industrial_demand.loc[spatial.ammonia.locations, "ammonia"].rename(index=lambda x: x + " NH3") / 8760
        else:
            p_set = industrial_demand["ammonia"].sum() / 8760

        n.madd("Load",
            spatial.ammonia.nodes,
            bus=spatial.ammonia.nodes,
            carrier="NH3",
            p_set=p_set
        )


def add_waste_heat(n,beccs):
    # TODO options?

    logger.info("Add possibility to use industrial waste heat in district heating")

    # AC buses with district heating
    urban_central = n.buses.index[n.buses.carrier == "urban central heat"]
    if not urban_central.empty:
        urban_central = urban_central.str[:-len(" urban central heat")]

        # TODO what is the 0.95 and should it be a config option?
        if options['use_fischer_tropsch_waste_heat']:
            n.links.loc[urban_central + " Fischer-Tropsch", "bus4"] = urban_central + " urban central heat"
            n.links.loc[urban_central + " Fischer-Tropsch", "efficiency4"] = (0.8 - n.links.loc[
                urban_central + " Fischer-Tropsch", "efficiency"]) * options['waste_heat_usage_share']

        for o in opts:
            if "B" in o:
                if options['use_biofuel_waste_heat']:
                    n.links.loc[urban_central + " biomass to liquid", "bus4"] = urban_central + " urban central heat"
                    n.links.loc[urban_central + " biomass to liquid", "efficiency4"] = (0.8 - n.links.loc[
                        urban_central + " biomass to liquid", "efficiency"]) * options['waste_heat_usage_share']

                    if beccs:
                        n.links.loc[urban_central + " biomass to liquid CC", "bus4"] = urban_central + " urban central heat"
                        n.links.loc[urban_central + " biomass to liquid CC", "efficiency4"] = (0.8 - n.links.loc[
                            urban_central + " biomass to liquid CC", "efficiency"]) * options['waste_heat_usage_share']

                    if options['electrobiofuels']:
                        n.links.loc[urban_central + " electrobiofuels", "bus4"] = urban_central + " urban central heat"
                        n.links.loc[urban_central + " electrobiofuels", "efficiency4"] = (0.8 - costs.at['electrobiofuels', "efficiency-tot"]) * options['waste_heat_usage_share']

                    n.links.loc[urban_central + " solid biomass to gas", "bus4"] = urban_central + " urban central heat"
                    n.links.loc[urban_central + " solid biomass to gas", "efficiency4"] = (0.8 - n.links.loc[
                        urban_central + " solid biomass to gas", "efficiency"]) * options['waste_heat_usage_share']

        if options['use_fuel_cell_waste_heat']:
            n.links.loc[urban_central + " H2 Fuel Cell", "bus4"] = urban_central + " urban central heat"
            n.links.loc[urban_central + " H2 Fuel Cell", "efficiency4"] = (0.8 - n.links.loc[
                urban_central + " H2 Fuel Cell", "efficiency"]) * options['waste_heat_usage_share']

        if options['use_electrolysis_waste_heat']:
            n.links.loc[urban_central + " H2 Electrolysis", "bus4"] = urban_central + " urban central heat"
            n.links.loc[urban_central + " H2 Electrolysis", "efficiency4"] = costs.at['electrolysis', 'efficiency-heat'] * options['waste_heat_usage_share']

        # heat_carriers = ["urban central heat", "services urban decentral heat"]
        # heat_buses = n.buses.index[n.buses.carrier.isin(heat_carriers)]
        # locations = n.buses.location[heat_buses]
        if options['use_dac_waste_heat']:
            n.links.loc[urban_central + " urban central DAC", "bus4"] = urban_central + " urban central heat"
            n.links.loc[urban_central + " urban central DAC", "efficiency4"] = costs.at['solid direct air capture', 'compression-heat-output'] * options['waste_heat_usage_share']
            # n.links.loc[urban_central + " services urban decentral DAC", "bus4"] = urban_central + " services urban decentral heat"
            # n.links.loc[urban_central + " services urban decentral DAC", "efficiency4"] = costs.at['direct air capture', 'compression-heat-output'] * options['waste_heat_usage_share']

def add_agriculture(n, costs):

    logger.info('Add agriculture, forestry and fishing sector.')

    nodes = pop_layout.index

    # electricity

    n.madd("Load",
        nodes,
        suffix=" agriculture electricity",
        bus=nodes,
        carrier='agriculture electricity',
        p_set=pop_weighted_energy_totals.loc[nodes, "total agriculture electricity"] * 1e6 / 8760
    )

    # heat

    n.madd("Load",
        nodes,
        suffix=" agriculture heat",
        bus=nodes + " services rural heat",
        carrier="agriculture heat",
        p_set=pop_weighted_energy_totals.loc[nodes, "total agriculture heat"] * 1e6 / 8760
    )

    # machinery

    electric_share = get(options["agriculture_machinery_electric_share"], investment_year)
    assert electric_share <= 1.
    ice_share = 1 - electric_share

    machinery_nodal_energy = pop_weighted_energy_totals.loc[nodes, "total agriculture machinery"]

    if electric_share > 0:

        efficiency_gain = options["agriculture_machinery_fuel_efficiency"] / options["agriculture_machinery_electric_efficiency"]

        n.madd("Load",
            nodes,
            suffix=" agriculture machinery electric",
            bus=nodes,
            carrier="agriculture machinery electric",
            p_set=electric_share / efficiency_gain * machinery_nodal_energy * 1e6 / 8760,
        )

    if ice_share > 0:

        n.madd("Load",
            ["agriculture machinery oil"],
            bus=spatial.oil.nodes,
            carrier="agriculture machinery oil",
            p_set=ice_share * machinery_nodal_energy.sum() * 1e6 / 8760
        )

        co2 = ice_share * machinery_nodal_energy.sum() * 1e6 / 8760 * costs.at["oil", 'CO2 intensity']

        n.add("Load",
            "agriculture machinery oil emissions",
            bus="co2 atmosphere",
            carrier="agriculture machinery oil emissions",
            p_set=-co2
        )

def decentral(n):
    """Removes the electricity transmission system."""
    n.lines.drop(n.lines.index, inplace=True)
    n.links.drop(n.links.index[n.links.carrier.isin(["DC", "B2B"])], inplace=True)


def remove_h2_network(n):

    n.links.drop(n.links.index[n.links.carrier.str.contains("H2 pipeline")], inplace=True)

    if "EU H2 Store" in n.stores.index:
        n.stores.drop("EU H2 Store", inplace=True)


def maybe_adjust_costs_and_potentials(n, opts):
    for o in opts:
        if "+" not in o: continue
        oo = o.split("+")
        carrier_list = np.hstack((n.generators.carrier.unique(), n.links.carrier.unique(),
                                  n.stores.carrier.unique(), n.storage_units.carrier.unique()))
        suptechs = map(lambda c: c.split("-", 2)[0], carrier_list)
        if oo[0].startswith(tuple(suptechs)):
            carrier = oo[0]
            attr_lookup = {"p": "p_nom_max", "e": "e_nom_max", "c": "capital_cost"}
            attr = attr_lookup[oo[1][0]]
            factor = float(oo[1][1:])
            # beware if factor is 0 and p_nom_max is np.inf, 0*np.inf is nan
            if carrier == "AC":  # lines do not have carrier
                n.lines[attr] *= factor
            else:
                if attr == 'p_nom_max':
                    comps = {"Generator", "Link", "StorageUnit"}
                elif attr == 'e_nom_max':
                    comps = {"Store"}
                else:
                    comps = {"Generator", "Link", "StorageUnit", "Store"}
                for c in n.iterate_components(comps):
                    if carrier == 'solar':
                        sel = c.df.carrier.str.contains(carrier) & ~c.df.carrier.str.contains("solar rooftop")
                    else:
                        sel = c.df.carrier.str.contains(carrier)
                    c.df.loc[sel, attr] *= factor
            print("changing", attr, "for", carrier, "by factor", factor)


def add_biomass_transport(n):
    # costs for biomass transport
    transport_costs = pd.read_csv('resources/biomass_transport_costs.csv',#snakemake.input.build_biomass_transport_costs,#biomass_transport,
                                  index_col=0)
    # add biomass transport
    biomass_transport = create_network_topology(n, "Biomass transport ")

    # make transport in both directions
    df = biomass_transport.copy()
    df["bus1"] = biomass_transport.bus0
    df["bus0"] = biomass_transport.bus1
    df.rename(index=lambda x: "Biomass transport " + df.at[x, "bus0"]
                              + " -> " + df.at[x, "bus1"], inplace=True)
    biomass_transport = pd.concat([biomass_transport, df])

    # costs
    bus0_costs = biomass_transport.bus0.apply(lambda x: transport_costs.loc[x[:2]])
    bus1_costs = biomass_transport.bus1.apply(lambda x: transport_costs.loc[x[:2]])
    biomass_transport["costs"] = pd.concat([bus0_costs, bus1_costs],
                                           axis=1).mean(axis=1)

    n.madd("Link",
           biomass_transport.index,
           bus0=biomass_transport.bus0 + " solid biomass",
           bus1=biomass_transport.bus1 + " solid biomass",
           p_nom_extendable=True,
           length=biomass_transport.length.values,
           marginal_cost=biomass_transport.costs * 0.01 * biomass_transport.length.values,
           capital_cost=1,
           carrier="solid biomass transport")


# TODO this should rather be a config no wildcard
def limit_individual_line_extension(n, maxext):
    logger.info(f"limiting new HVAC and HVDC extensions to {maxext} MW")
    n.lines['s_nom_max'] = n.lines['s_nom'] + maxext
    hvdc = n.links.index[n.links.carrier == 'DC']
    n.links.loc[hvdc, 'p_nom_max'] = n.links.loc[hvdc, 'p_nom'] + maxext


def hvdc_transport_model(n):
    print("Changing AC lines to HVDC links")
    n.madd("Link",
           n.lines.index,
           bus0=n.lines.bus0,
           bus1=n.lines.bus1,
           p_nom_extendable=True,
           p_nom=n.lines.s_nom,
           p_nom_min=n.lines.s_nom,
           p_min_pu=-1,
           efficiency=1 - 0.03 * n.lines.length / 1000,
           marginal_cost=0,
           carrier="DC",
           length=n.lines.length,
           capital_cost=n.lines.capital_cost)

    # Remove AC lines
    print("Removing AC lines")
    lines_rm = n.lines.index
    n.mremove("Line", lines_rm)

    # Set efficiency of all DC links to include losses depending on length
    n.links.loc[n.links.carrier == 'DC', 'efficiency'] = 1 - 0.03 * n.links.loc[
        n.links.carrier == 'DC', 'length'] / 1000

def no_transmission(n):
    # Remove AC lines
    print("Removing AC lines")
    lines_rm = n.lines.index
    n.mremove("Line", lines_rm)

    n.links.loc[n.links.carrier == 'DC', 'p_nom'] = 0
    # n.links.loc[n.links.carrier == 'DC', 'p_nom_extendable'] = 'False'

def apply_time_segmentation(n, segments, solver_name="cbc",
                            overwrite_time_dependent=True):
    """Aggregating time series to segments with different lengths

    Input:
        n: pypsa Network
        segments: (int) number of segments in which the typical period should be
                  subdivided
        solver_name: (str) name of solver
        overwrite_time_dependent: (bool) overwrite time dependent data of pypsa network
        with typical time series created by tsam
    """
    try:
        import tsam.timeseriesaggregation as tsam
    except:
        raise ModuleNotFoundError("Optional dependency 'tsam' not found."
                                  "Install via 'pip install tsam'")

    # get all time-dependent data
    columns = pd.MultiIndex.from_tuples([],names=['component', 'key', 'asset'])
    raw = pd.DataFrame(index=n.snapshots,columns=columns)
    for c in n.iterate_components():
        for attr, pnl in c.pnl.items():
            # exclude e_min_pu which is used for SOC of EVs in the morning
            if not pnl.empty and attr != 'e_min_pu':
                df = pnl.copy()
                df.columns = pd.MultiIndex.from_product([[c.name], [attr], df.columns])
                raw = pd.concat([raw, df], axis=1)

    # normalise all time-dependent data
    annual_max = raw.max().replace(0,1)
    raw = raw.div(annual_max, level=0)

    # get representative segments
    agg = tsam.TimeSeriesAggregation(raw, hoursPerPeriod=len(raw),
                                     noTypicalPeriods=1, noSegments=int(segments),
                                     segmentation=True, solver=solver_name)
    segmented = agg.createTypicalPeriods()


    weightings = segmented.index.get_level_values("Segment Duration")
    offsets = np.insert(np.cumsum(weightings[:-1]), 0, 0)
    timesteps = [raw.index[0] + pd.Timedelta(f"{offset}h") for offset in offsets]
    snapshots =  pd.DatetimeIndex(timesteps)
    sn_weightings = pd.Series(weightings, index=snapshots, name="weightings", dtype="float64")

    n.set_snapshots(sn_weightings.index)
    n.snapshot_weightings = n.snapshot_weightings.mul(sn_weightings, axis=0)

    # overwrite time-dependent data with timeseries created by tsam
    if overwrite_time_dependent:
        values_t = segmented.mul(annual_max).set_index(snapshots)
        for component, key in values_t.columns.droplevel(2).unique():
            n.pnl(component)[key] = values_t[component, key]

    return n

def set_temporal_aggregation(n, opts, solver_name):
    """Aggregate network temporally."""
    for o in opts:
        # temporal averaging
        m = re.match(r"^\d+h$", o, re.IGNORECASE)
        if m is not None:
            n = average_every_nhours(n, m.group(0))
            break
        # representive snapshots
        m = re.match(r"(^\d+)sn$", o, re.IGNORECASE)
        if m is not None:
            sn = int(m[1])
            logger.info(f"use every {sn} snapshot as representative")
            n.set_snapshots(n.snapshots[::sn])
            n.snapshot_weightings *= sn
            break
        # segments with package tsam
        m = re.match(r"^(\d+)seg$", o, re.IGNORECASE)
        if m is not None:
            segments = int(m[1])
            logger.info(f"use temporal segmentation with {segments} segments")
            n = apply_time_segmentation(n, segments, solver_name=solver_name)
            break
    return n

#%%

if __name__ == "__main__":
    if 'snakemake' not in globals():
        from helper import mock_snakemake

        snakemake = mock_snakemake(
            'prepare_sector_network',
            simpl='',
            opts="",
            clusters="37",
            lv=1.5,
            sector_opts='cb40ex0-365H-T-H-B-I-A-solar+p3-dist1',
            planning_horizons="2020",
        )

    logging.basicConfig(level=snakemake.config['logging_level'])

    update_config_with_sector_opts(snakemake.config, snakemake.wildcards.sector_opts)

    options = snakemake.config["sector"]

    opts = snakemake.wildcards.sector_opts.split('-')

    print('Options: ', opts)

    investment_year = int(snakemake.wildcards.planning_horizons[-4:])

    overrides = override_component_attrs(snakemake.input.overrides)
    n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    Nyears = n.snapshot_weightings.generators.sum() / 8760

    costs = prepare_costs(snakemake.input.costs,
                          snakemake.config['costs']['USD2013_to_EUR2013'],
                          snakemake.config['costs']['discountrate'],
                          Nyears,
                          snakemake.config['costs']['lifetime'])

    biomass_import_price = snakemake.config['biomass']['import price'] * 3.6  # EUR/MWh
    carbon_sequestration_cost = options["co2_sequestration_cost"]
    costs, biomass_import_price, carbon_sequestration_cost, options = sensitivity_costs(costs, biomass_import_price, carbon_sequestration_cost, options)

    print('Updated biomass import price: ', biomass_import_price)
    print('Updated CO2 sequestration cost: ', carbon_sequestration_cost)

    pop_weighted_energy_totals = pd.read_csv(snakemake.input.pop_weighted_energy_totals, index_col=0)

    patch_electricity_network(n)

    spatial = define_spatial(pop_layout.index, options)

    if snakemake.config["foresight"] == 'myopic':

        add_lifetime_wind_solar(n, costs)

        conventional = snakemake.config['existing_capacities']['conventional_carriers']
        for carrier in conventional:
            add_carrier_buses(n, carrier)

    add_co2_tracking(n, options, carbon_sequestration_cost)

    add_generation(n, costs)

    add_storage_and_grids(n, costs)

    # TODO merge with opts cost adjustment below
    for o in opts:
        if o[:4] == "wave":
            wave_cost_factor = float(o[4:].replace("p", ".").replace("m", "-"))
            print("Including wave generators with cost factor of", wave_cost_factor)
            add_wave(n, wave_cost_factor)
        if o[:4] == "dist":
            options['electricity_distribution_grid'] = True
            options['electricity_distribution_grid_cost_factor'] = float(o[4:].replace("p", ".").replace("m", "-"))

    if "nodistrict" in opts:
        options["district_heating"]["progress"] = 0.0

    if "T" in opts:
        add_land_transport(n, costs)

    if "H" in opts:
        add_heat(n, costs)

    if "I" in opts:
        add_industry(n, costs)

    if options["hvdc"] and not options["notransmission"]:
        hvdc_transport_model(n)

    if options["notransmission"]:
        no_transmission(n)

    beccs = snakemake.config['biomass']['beccs']
    for o in opts:
        if "B" in o:
            add_biomass(n, costs, beccs, biomass_import_price)
            if options["biomass_transport"]:
                add_biomass_transport(n)

    if options['dac']:
        add_dac(n, costs)

    if "I" in opts and "H" in opts:
        add_waste_heat(n,beccs)

    if "A" in opts:  # requires H and I
        add_agriculture(n, costs)


    if options['ammonia']:
        add_ammonia(n, costs)

    if "decentral" in opts:
        decentral(n)

    if "noH2network" in opts:
        remove_h2_network(n)

    if options["co2_network"]:
        add_co2_network(n, costs)

    solver_name = snakemake.config["solving"]["solver"]["name"]
    n = set_temporal_aggregation(n, opts, solver_name)

    limit_type = "config"
    limit = get(snakemake.config["co2_budget"], investment_year)
    for o in opts:
        if not "cb" in o: continue
        limit_type = "carbon budget"
        fn = snakemake.config['results_dir'] + snakemake.config['run'] + '/csvs/carbon_budget_distribution.csv'
        if not os.path.exists(fn):
            emissions_scope = snakemake.config["energy"]["emissions"]
            report_year = snakemake.config["energy"]["eurostat_report_year"]
            build_carbon_budget(o, snakemake.input.eurostat, fn, emissions_scope, report_year)
        co2_cap = pd.read_csv(fn, index_col=0).squeeze()
        limit = co2_cap.loc[investment_year]
        break
    for o in opts:
        if not "Co2L" in o: continue
        limit_type = "wildcard"
        limit = o[o.find("Co2L") + 4:]
        limit = float(limit.replace("p", ".").replace("m", "-"))
        break
    print("Add CO2 limit from", limit_type)
    add_co2limit(n, Nyears, limit)

    for o in opts:
        if not o[:10] == 'linemaxext': continue
        maxext = float(o[10:]) * 1e3
        limit_individual_line_extension(n, maxext)
        break

    if options['electricity_distribution_grid']:
        insert_electricity_distribution_grid(n, costs)

    maybe_adjust_costs_and_potentials(n, opts)

    if options['gas_distribution_grid']:
        insert_gas_distribution_costs(n, costs)

    if options['electricity_grid_connection']:
        add_electricity_grid_connection(n, costs)
        
    if options["Ethylene"]:
    	add_ethylene(n,costs)
    
    if options["EO"]:
        add_EO(n)
    	
    if options["MEA"]:
    	add_MEA(n)
    	
    if options["KOH"]:
        add_KOH(n)
        
    if options["PEI silica sorbent"]:
        add_PEI(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output[0])
