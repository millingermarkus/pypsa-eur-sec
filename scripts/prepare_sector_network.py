# coding: utf-8

import logging
logger = logging.getLogger(__name__)
import pandas as pd
idx = pd.IndexSlice

import numpy as np
import xarray as xr
import re, os, sys

from six import iteritems, string_types

import pypsa

import yaml

import pytz

from vresutils.costdata import annuity

from scipy.stats import beta
from build_energy_totals import build_eea_co2, build_eurostat_co2, build_co2_totals

#First tell PyPSA that links can have multiple outputs by
#overriding the component_attrs. This can be done for
#as many buses as you need with format busi for i = 2,3,4,5,....
#See https://pypsa.org/doc/components.html#link-with-multiple-outputs-or-inputs


override_component_attrs = pypsa.descriptors.Dict({k : v.copy() for k,v in pypsa.components.component_attrs.items()})
override_component_attrs["Link"].loc["bus2"] = ["string",np.nan,np.nan,"2nd bus","Input (optional)"]
override_component_attrs["Link"].loc["bus3"] = ["string",np.nan,np.nan,"3rd bus","Input (optional)"]
override_component_attrs["Link"].loc["bus4"] = ["string",np.nan,np.nan,"4th bus","Input (optional)"]
override_component_attrs["Link"].loc["efficiency2"] = ["static or series","per unit",1.,"2nd bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["efficiency3"] = ["static or series","per unit",1.,"3rd bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["efficiency4"] = ["static or series","per unit",1.,"4th bus efficiency","Input (optional)"]
override_component_attrs["Link"].loc["p2"] = ["series","MW",0.,"2nd bus output","Output"]
override_component_attrs["Link"].loc["p3"] = ["series","MW",0.,"3rd bus output","Output"]
override_component_attrs["Link"].loc["p4"] = ["series","MW",0.,"4th bus output","Output"]

override_component_attrs["Link"].loc["build_year"] = ["integer","year",np.nan,"build year","Input (optional)"]
override_component_attrs["Link"].loc["lifetime"] = ["float","years",np.nan,"lifetime","Input (optional)"]
override_component_attrs["Generator"].loc["build_year"] = ["integer","year",np.nan,"build year","Input (optional)"]
override_component_attrs["Generator"].loc["lifetime"] = ["float","years",np.nan,"lifetime","Input (optional)"]
override_component_attrs["Store"].loc["build_year"] = ["integer","year",np.nan,"build year","Input (optional)"]
override_component_attrs["Store"].loc["lifetime"] = ["float","years",np.nan,"lifetime","Input (optional)"]



def co2_emissions_year(cts, opts, year):
    """
    Calculate CO2 emissions in one specific year (e.g. 1990 or 2018).
    """
    eea_co2 = build_eea_co2(year)

    # TODO: read Eurostat data from year>2014, this only affects the estimation of
    # CO2 emissions for "BA","RS","AL","ME","MK"
    if year > 2014:
        eurostat_co2 = build_eurostat_co2(year=2014)
    else:
        eurostat_co2 = build_eurostat_co2(year)

    co2_totals = build_co2_totals(eea_co2, eurostat_co2)

    co2_emissions = co2_totals.loc[cts, "electricity"].sum()

    if "T" in opts:
        co2_emissions += co2_totals.loc[cts, [i+ " non-elec" for i in ["rail","road"]]].sum().sum()
    if "H" in opts:
        co2_emissions += co2_totals.loc[cts, [i+ " non-elec" for i in ["residential","services"]]].sum().sum()
    if "I" in opts:
        co2_emissions += co2_totals.loc[cts, ["industrial non-elec","industrial processes",
                                              "domestic aviation","international aviation",
                                              "domestic navigation","international navigation"]].sum().sum()
	
    co2_emissions *= 0.001  # Convert MtCO2 to GtCO2
    return co2_emissions


def build_carbon_budget(o):
    #distribute carbon budget following beta or exponential transition path
    if "be" in o:
        #beta decay
        carbon_budget = float(o[o.find("cb")+2:o.find("be")])
        be=float(o[o.find("be")+2:])
    if "ex" in o:
        #exponential decay
        carbon_budget = float(o[o.find("cb")+2:o.find("ex")])
        r=float(o[o.find("ex")+2:])


    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)
    pop_layout["ct"] = pop_layout.index.str[:2]
    cts = pop_layout.ct.value_counts().index

    e_1990 = co2_emissions_year(cts, opts, year=1990)

    #emissions at the beginning of the path (last year available 2018)
    e_0 = co2_emissions_year(cts, opts, year=2018)
    #emissions in 2019 and 2020 assumed equal to 2018 and substracted
    carbon_budget -= 2*e_0
    planning_horizons = snakemake.config['scenario']['planning_horizons']
    CO2_CAP = pd.DataFrame(index = pd.Series(data=planning_horizons,
                                             name='planning_horizon'),
                                             columns=pd.Series(data=[],
                                                              name='paths',
                                                              dtype='float'))
    t_0 = planning_horizons[0]
    if "be" in o:
        #beta decay
        t_f = t_0 + (2*carbon_budget/e_0).round(0) # final year in the path
        #emissions (relative to 1990)
        CO2_CAP[o] = [(e_0/e_1990)*(1-beta.cdf((t-t_0)/(t_f-t_0), be, be)) for t in planning_horizons]

    if "ex" in o:
        #exponential decay without delay
        T=carbon_budget/e_0
        m=(1+np.sqrt(1+r*T))/T
        CO2_CAP[o] = [(e_0/e_1990)*(1+(m+r)*(t-t_0))*np.exp(-m*(t-t_0)) for t in planning_horizons]


    CO2_CAP.to_csv(path_cb + 'carbon_budget_distribution.csv', sep=',',
                   line_terminator='\n',  float_format='%.3f')
    countries=pd.Series(data=cts)
    countries.to_csv(path_cb + 'countries.csv', sep=',',
               line_terminator='\n',  float_format='%.3f')


def add_lifetime_wind_solar(n):
    """
    Add lifetime for solar and wind generators
    """
    for carrier in ['solar', 'onwind', 'offwind-dc', 'offwind-ac']:
        carrier_name='offwind' if carrier in ['offwind-dc', 'offwind-ac'] else carrier
        n.generators.loc[[index for index in n.generators.index.to_list()
                        if carrier in index], 'lifetime']=costs.at[carrier_name,'lifetime']


def create_network_topology(n, prefix="", connector=" -> "):
    """
    create a network topology as the electric network,
    returns a pandas dataframe with bus0, bus1 and length
    """
    attrs = ["bus0", "bus1", "length"]

    candidates = pd.concat([n.lines[attrs],
                            n.links.loc[n.links.carrier == "DC", attrs]])

    positive_order = candidates.bus0 < candidates.bus1
    candidates_p = candidates[positive_order]
    candidates_n = (candidates[~ positive_order]
                    .rename(columns={"bus0": "bus1", "bus1": "bus0"}))
    candidates = pd.concat((candidates_p, candidates_n), sort=False)

    topo = candidates.groupby(["bus0", "bus1"], as_index=False).mean()
    topo.index = topo.apply(lambda x: prefix + x.bus0 + connector + x.bus1, axis=1)
    return topo


def update_wind_solar_costs(n,costs):
    """
    Update costs for wind and solar generators added with pypsa-eur to those
    cost in the planning year

    """

    #NB: solar costs are also manipulated for rooftop
    #when distribution grid is inserted
    n.generators.loc[n.generators.carrier=='solar','capital_cost'] = costs.at['solar-utility', 'fixed']

    n.generators.loc[n.generators.carrier=='onwind','capital_cost'] = costs.at['onwind', 'fixed']

    #for offshore wind, need to calculated connection costs

    #assign clustered bus
    #map initial network -> simplified network
    busmap_s = pd.read_csv(snakemake.input.busmap_s, index_col=0).squeeze()
    busmap_s.index = busmap_s.index.astype(str)
    busmap_s = busmap_s.astype(str)
    #map simplified network -> clustered network
    busmap = pd.read_csv(snakemake.input.busmap, index_col=0).squeeze()
    busmap.index = busmap.index.astype(str)
    busmap = busmap.astype(str)
    #map initial network -> clustered network
    clustermaps = busmap_s.map(busmap)

    #code adapted from pypsa-eur/scripts/add_electricity.py
    for connection in ['dc','ac']:
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

            #convert to aggregated clusters with weighting
            weight = ds['weight'].to_pandas()

            #e.g. clusters == 37m means that VRE generators are left
            #at clustering of simplified network, but that they are
            #connected to 37-node network
            if snakemake.wildcards.clusters[-1:] == "m":
                genmap = busmap_s
            else:
                genmap = clustermaps

            connection_cost = (connection_cost*weight).groupby(genmap).sum()/weight.groupby(genmap).sum()

            capital_cost = (costs.at['offwind', 'fixed'] +
                            costs.at[tech + '-station', 'fixed'] +
                            connection_cost)

            logger.info("Added connection cost of {:0.0f}-{:0.0f} Eur/MW/a to {}"
                        .format(connection_cost[0].min(), connection_cost[0].max(), tech))

            n.generators.loc[n.generators.carrier==tech,'capital_cost'] = capital_cost.rename(index=lambda node: node + ' ' + tech)


def add_carrier_buses(n, carriers):
    """
    Add buses to connect e.g. coal, nuclear and oil plants
    """

    for carrier in carriers:

        n.add("Carrier",
              carrier)

        #use madd to get location inserted
        n.madd("Bus",
               ["EU " + carrier],
               location="EU",
               carrier=carrier)

        #use madd to get carrier inserted
        n.madd("Store",
               ["EU " + carrier + " Store"],
               bus=["EU " + carrier],
               e_nom_extendable=True,
               e_cyclic=True,
               carrier=carrier,
               capital_cost=0.) #could correct to e.g. 0.2 EUR/kWh * annuity and O&M

        n.add("Generator",
              "EU " + carrier,
              bus="EU " + carrier,
              p_nom_extendable=True,
              carrier=carrier,
              capital_cost=0.,
              marginal_cost=costs.at[carrier,'fuel'])


def remove_elec_base_techs(n):
    """remove conventional generators (e.g. OCGT) and storage units (e.g. batteries and H2)
    from base electricity-only network, since they're added here differently using links
    """

    for c in n.iterate_components(snakemake.config["pypsa_eur"]):
        to_keep = snakemake.config["pypsa_eur"][c.name]
        to_remove = pd.Index(c.df.carrier.unique()).symmetric_difference(to_keep)
        print("Removing",c.list_name,"with carrier",to_remove)
        names = c.df.index[c.df.carrier.isin(to_remove)]
        print(names)
        n.mremove(c.name, names)
        n.carriers.drop(to_remove, inplace=True, errors="ignore")


def remove_non_electric_buses(n):
    """
    remove buses from pypsa-eur with carriers which are not AC buses
    """
    print("drop buses from PyPSA-Eur with carrier: ", n.buses[~n.buses.carrier.isin(["AC", "DC"])].carrier.unique())
    n.buses = n.buses[n.buses.carrier.isin(["AC", "DC"])]


def add_co2_tracking(n):


    #minus sign because opposite to how fossil fuels used:
    #CH4 burning puts CH4 down, atmosphere up
    n.add("Carrier","co2",
          co2_emissions=-1.)

    #this tracks CO2 in the atmosphere
    n.madd("Bus",
           ["co2 atmosphere"],
           location="EU",
           carrier="co2")

    #NB: can also be negative
    n.madd("Store",["co2 atmosphere"],
           e_nom_extendable=True,
           e_min_pu=-1,
           carrier="co2",
           bus="co2 atmosphere")

    #this tracks CO2 stored, e.g. underground
    n.madd("Bus",
           ["co2 stored"],
           location="EU",
           carrier="co2 stored")

    co2_sequestration_potential = 0
    for o in opts:
        if "S" in o:
            co2_sequestration_potential = float(o[o.find("S") + 1:])

    print('CO2 sequestration potential: ', co2_sequestration_potential)

    n.madd("Store",["co2 stored"],
           e_nom_extendable=True,
           e_nom_max=co2_sequestration_potential*1e6, #options['co2_sequestration_potential']*1e6,
           capital_cost=options['co2_sequestration_cost'],
           carrier="co2 stored",
           bus="co2 stored")

    if options['co2_vent']:
        n.madd("Link",["co2 vent"],
               bus0="co2 stored",
               bus1="co2 atmosphere",
               carrier="co2 vent",
               efficiency=1.,
               p_nom_extendable=True)

def add_dac(n):

    heat_buses = n.buses.index[n.buses.carrier.isin(["urban central heat",
                                                     "services urban decentral heat"])]
    locations = n.buses.location[heat_buses]

    n.madd("Link",
           locations,
           suffix=" DAC",
           bus0="co2 atmosphere",
           bus1="co2 stored",
           bus2=locations.values,
           bus3=heat_buses,
           carrier="DAC",
           capital_cost=costs.at['direct air capture','fixed'],
           efficiency=1.,
           efficiency2=-(costs.at['direct air capture','electricity-input'] + costs.at['direct air capture','compression-electricity-input']),
           efficiency3=-(costs.at['direct air capture','heat-input'] - costs.at['direct air capture','compression-heat-output']),
           p_nom_extendable=True,
           lifetime=costs.at['direct air capture','lifetime'])


def add_co2limit(n, Nyears=1.,limit=0.):

    cts = pop_layout.ct.value_counts().index

    co2_limit = co2_totals.loc[cts, "electricity"].sum()

    if "T" in opts:
        co2_limit += co2_totals.loc[cts, [i+ " non-elec" for i in ["rail","road"]]].sum().sum()
    if "H" in opts:
        co2_limit += co2_totals.loc[cts, [i+ " non-elec" for i in ["residential","services"]]].sum().sum()
    if "I" in opts:
        co2_limit += co2_totals.loc[cts, ["industrial non-elec","industrial processes",
                                          "domestic aviation","international aviation",
                                          "domestic navigation","international navigation"]].sum().sum()

    co2_limit *= limit*Nyears

    n.add("GlobalConstraint", "CO2Limit",
          carrier_attribute="co2_emissions", sense="<=",
          constant=co2_limit)

def add_emission_prices(n, emission_prices=None, exclude_co2=False):
    assert False, "Needs to be fixed, adds NAN"

    if emission_prices is None:
        emission_prices = snakemake.config['costs']['emission_prices']
    if exclude_co2: emission_prices.pop('co2')
    ep = (pd.Series(emission_prices).rename(lambda x: x+'_emissions') * n.carriers).sum(axis=1)
    n.generators['marginal_cost'] += n.generators.carrier.map(ep)
    n.storage_units['marginal_cost'] += n.storage_units.carrier.map(ep)

def set_line_s_max_pu(n):
    # set n-1 security margin to 0.5 for 37 clusters and to 0.7 from 200 clusters
    # 128 reproduces 98% of line volume in TWkm, but clustering distortions inside node
    n_clusters = len(n.buses.index[n.buses.carrier == "AC"])
    s_max_pu = np.clip(0.5 + 0.2 * (n_clusters - 37) / (200 - 37), 0.5, 0.7)
    n.lines['s_max_pu'] = s_max_pu

    dc_b = n.links.carrier == 'DC'
    n.links.loc[dc_b, 'p_max_pu'] = snakemake.config['links']['p_max_pu']
    n.links.loc[dc_b, 'p_min_pu'] = - snakemake.config['links']['p_max_pu']

def set_line_volume_limit(n, lv):

    dc_b = n.links.carrier == 'DC'

    if lv != "opt":
        lv = float(lv)

        # Either line_volume cap or cost
        n.lines['capital_cost'] = 0.
        n.links.loc[dc_b,'capital_cost'] = 0.
    else:
        n.lines['capital_cost'] = (n.lines['length'] *
                                   costs.at['HVAC overhead', 'fixed'])

        #add HVDC inverter post factor, to maintain consistency with LV limit
        n.links.loc[dc_b, 'capital_cost'] = (n.links.loc[dc_b, 'length'] *
                                             costs.at['HVDC overhead', 'fixed'])# +
                                             #costs.at['HVDC inverter pair', 'fixed'])



    if lv != 1.0:
        lines_s_nom = n.lines.s_nom.where(
            n.lines.type == '',
            np.sqrt(3) * n.lines.num_parallel *
            n.lines.type.map(n.line_types.i_nom) *
            n.lines.bus0.map(n.buses.v_nom)
        )

        n.lines['s_nom_min'] = lines_s_nom

        n.links.loc[dc_b,'p_nom_min'] = n.links['p_nom']

        n.lines['s_nom_extendable'] = True
        n.links.loc[dc_b,'p_nom_extendable'] = True

        if lv != "opt":
            n.line_volume_limit = lv * ((lines_s_nom * n.lines['length']).sum() +
                                        n.links.loc[dc_b].eval('p_nom * length').sum())

    return n



# def _add_links_from_tyndp(buses, links):
#     links_tyndp = pd.read_csv(snakemake.input.links_tyndp)
#
#     # # remove all links from list which lie outside all of the desired countries
#     # europe_shape = gpd.read_file(snakemake.input.europe_shape).loc[0, 'geometry']
#     # europe_shape_prepped = shapely.prepared.prep(europe_shape)
#     # x1y1_in_europe_b = links_tyndp[['x1', 'y1']].apply(lambda p: europe_shape_prepped.contains(Point(p)), axis=1)
#     # x2y2_in_europe_b = links_tyndp[['x2', 'y2']].apply(lambda p: europe_shape_prepped.contains(Point(p)), axis=1)
#     # is_within_covered_countries_b = x1y1_in_europe_b & x2y2_in_europe_b
#     #
#     # if not is_within_covered_countries_b.all():
#     #     logger.info("TYNDP links outside of the covered area (skipping): " +
#     #                 ", ".join(links_tyndp.loc[~ is_within_covered_countries_b, "Name"]))
#     #
#     #     links_tyndp = links_tyndp.loc[is_within_covered_countries_b]
#     #     if links_tyndp.empty:
#     #         return buses, links
#
#     has_replaces_b = links_tyndp.replaces.notnull()
#     oids = dict(Bus=_get_oid(buses), Link=_get_oid(links))
#     keep_b = dict(Bus=pd.Series(True, index=buses.index),
#                   Link=pd.Series(True, index=links.index))
#     for reps in links_tyndp.loc[has_replaces_b, 'replaces']:
#         for comps in reps.split(':'):
#             oids_to_remove = comps.split('.')
#             c = oids_to_remove.pop(0)
#             keep_b[c] &= ~oids[c].isin(oids_to_remove)
#     buses = buses.loc[keep_b['Bus']]
#     links = links.loc[keep_b['Link']]
#
#     links_tyndp["j"] = _find_closest_links(links, links_tyndp, distance_upper_bound=0.20)
#     # Corresponds approximately to 20km tolerances
#
#     if links_tyndp["j"].notnull().any():
#         logger.info("TYNDP links already in the dataset (skipping): " + ", ".join(links_tyndp.loc[links_tyndp["j"].notnull(), "Name"]))
#         links_tyndp = links_tyndp.loc[links_tyndp["j"].isnull()]
#
#     tree = sp.spatial.KDTree(buses[['x', 'y']])
#     _, ind0 = tree.query(links_tyndp[["x1", "y1"]])
#     ind0_b = ind0 < len(buses)
#     links_tyndp.loc[ind0_b, "bus0"] = buses.index[ind0[ind0_b]]
#
#     _, ind1 = tree.query(links_tyndp[["x2", "y2"]])
#     ind1_b = ind1 < len(buses)
#     links_tyndp.loc[ind1_b, "bus1"] = buses.index[ind1[ind1_b]]
#
#     links_tyndp_located_b = links_tyndp["bus0"].notnull() & links_tyndp["bus1"].notnull()
#     if not links_tyndp_located_b.all():
#         logger.warning("Did not find connected buses for TYNDP links (skipping): " + ", ".join(links_tyndp.loc[~links_tyndp_located_b, "Name"]))
#         links_tyndp = links_tyndp.loc[links_tyndp_located_b]
#
#     logger.info("Adding the following TYNDP links: " + ", ".join(links_tyndp["Name"]))
#
#     links_tyndp = links_tyndp[["bus0", "bus1"]].assign(
#         carrier='DC',
#         p_nom=links_tyndp["Power (MW)"],
#         length=links_tyndp["Length (given) (km)"].fillna(links_tyndp["Length (distance*1.2) (km)"]),
#         under_construction=True,
#         underground=False,
#         geometry=(links_tyndp[["x1", "y1", "x2", "y2"]]
#                   .apply(lambda s: str(LineString([[s.x1, s.y1], [s.x2, s.y2]])), axis=1)),
#         tags=('"name"=>"' + links_tyndp["Name"] + '", ' +
#               '"ref"=>"' + links_tyndp["Ref"] + '", ' +
#               '"status"=>"' + links_tyndp["status"] + '"')
#     )
#
#     links_tyndp.index = "T" + links_tyndp.index.astype(str)
#
#     return buses, links.append(links_tyndp, sort=True)



def average_every_nhours(n, offset):
    logger.info('Resampling the network to {}'.format(offset))
    m = n.copy(with_time=False)

    #fix copying of network attributes
    #copied from pypsa/io.py, should be in pypsa/components.py#Network.copy()
    allowed_types = (float,int,bool) + string_types + tuple(np.typeDict.values())
    attrs = dict((attr, getattr(n, attr))
                 for attr in dir(n)
                 if (not attr.startswith("__") and
                     isinstance(getattr(n,attr), allowed_types)))
    for k,v in iteritems(attrs):
        setattr(m,k,v)

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name+"_t")
        for k, df in iteritems(c.pnl):
            if not df.empty:
                if c.list_name == "stores" and k == "e_max_pu":
                    pnl[k] = df.resample(offset).min()
                elif c.list_name == "stores" and k == "e_min_pu":
                    pnl[k] = df.resample(offset).max()
                else:
                    pnl[k] = df.resample(offset).mean()

    return m


def generate_periodic_profiles(dt_index=pd.date_range("2011-01-01 00:00","2011-12-31 23:00",freq="H",tz="UTC"),
                               nodes=[],
                               weekly_profile=range(24*7)):
    """Give a 24*7 long list of weekly hourly profiles, generate this for
       each country for the period dt_index, taking account of time
       zones and Summer Time.

    """


    weekly_profile = pd.Series(weekly_profile,range(24*7))

    week_df = pd.DataFrame(index=dt_index,columns=nodes)

    for ct in nodes:
        week_df[ct] = [24*dt.weekday()+dt.hour for dt in dt_index.tz_convert(pytz.timezone(timezone_mappings[ct[:2]]))]
        week_df[ct] = week_df[ct].map(weekly_profile)

    return week_df



def shift_df(df,hours=1):
    """Works both on Series and DataFrame"""
    df = df.copy()
    df.values[:] = np.concatenate([df.values[-hours:],
                                   df.values[:-hours]])
    return df

def transport_degree_factor(temperature,deadband_lower=15,deadband_upper=20,
                            lower_degree_factor=0.5,
                            upper_degree_factor=1.6):

    """Work out how much energy demand in vehicles increases due to heating and cooling.

    There is a deadband where there is no increase.

    Degree factors are % increase in demand compared to no heating/cooling fuel consumption.

    Returns per unit increase in demand for each place and time
    """

    dd = temperature.copy()

    dd[(temperature > deadband_lower) & (temperature < deadband_upper)] = 0.

    dd[temperature < deadband_lower] = lower_degree_factor/100.*(deadband_lower-temperature[temperature < deadband_lower])

    dd[temperature > deadband_upper] = upper_degree_factor/100.*(temperature[temperature > deadband_upper]-deadband_upper)

    return dd


def prepare_data(network):


    ##############
    #Heating
    ##############


    ashp_cop = xr.open_dataarray(snakemake.input.cop_air_total).T.to_pandas().reindex(index=network.snapshots)
    gshp_cop = xr.open_dataarray(snakemake.input.cop_soil_total).T.to_pandas().reindex(index=network.snapshots)

    solar_thermal = xr.open_dataarray(snakemake.input.solar_thermal_total).T.to_pandas().reindex(index=network.snapshots)
    #1e3 converts from W/m^2 to MW/(1000m^2) = kW/m^2
    solar_thermal = options['solar_cf_correction'] * solar_thermal/1e3

    energy_totals = pd.read_csv(snakemake.input.energy_totals_name,index_col=0)

    nodal_energy_totals = energy_totals.loc[pop_layout.ct].fillna(0.)
    nodal_energy_totals.index = pop_layout.index
    nodal_energy_totals = nodal_energy_totals.multiply(pop_layout.fraction,axis=0)

    #copy forward the daily average heat demand into each hour, so it can be multipled by the intraday profile
    daily_space_heat_demand = xr.open_dataarray(snakemake.input.heat_demand_total).T.to_pandas().reindex(index=network.snapshots, method="ffill")

    intraday_profiles = pd.read_csv(snakemake.input.heat_profile,index_col=0)

    sectors = ["residential","services"]
    uses = ["water","space"]

    heat_demand = {}
    electric_heat_supply = {}
    for sector in sectors:
        for use in uses:
            intraday_year_profile = generate_periodic_profiles(daily_space_heat_demand.index.tz_localize("UTC"),
                                                               nodes=daily_space_heat_demand.columns,
                                                               weekly_profile=(list(intraday_profiles["{} {} weekday".format(sector,use)])*5 + list(intraday_profiles["{} {} weekend".format(sector,use)])*2)).tz_localize(None)

            if use == "space":
                heat_demand_shape = daily_space_heat_demand*intraday_year_profile
            else:
                heat_demand_shape = intraday_year_profile

            heat_demand["{} {}".format(sector,use)] = (heat_demand_shape/heat_demand_shape.sum()).multiply(nodal_energy_totals["total {} {}".format(sector,use)])*1e6
            electric_heat_supply["{} {}".format(sector,use)] = (heat_demand_shape/heat_demand_shape.sum()).multiply(nodal_energy_totals["electricity {} {}".format(sector,use)])*1e6

    heat_demand = pd.concat(heat_demand,axis=1)
    electric_heat_supply = pd.concat(electric_heat_supply,axis=1)

    #subtract from electricity load since heat demand already in heat_demand
    electric_nodes = n.loads.index[n.loads.carrier == "electricity"]
    n.loads_t.p_set[electric_nodes] = n.loads_t.p_set[electric_nodes] - electric_heat_supply.groupby(level=1,axis=1).sum()[electric_nodes]

    ##############
    #Transport
    ##############


    ## Get overall demand curve for all vehicles

    traffic = pd.read_csv(os.path.join(snakemake.input.traffic_data, "KFZ__count"),
                          skiprows=2)["count"]

    #Generate profiles
    transport_shape = generate_periodic_profiles(dt_index=network.snapshots.tz_localize("UTC"),
                                                 nodes=pop_layout.index,
                                                 weekly_profile=traffic.values).tz_localize(None)
    transport_shape = transport_shape/transport_shape.sum()

    transport_data = pd.read_csv(snakemake.input.transport_name,
                                 index_col=0)

    nodal_transport_data = transport_data.loc[pop_layout.ct].fillna(0.)
    nodal_transport_data.index = pop_layout.index
    nodal_transport_data["number cars"] = pop_layout["fraction"]*nodal_transport_data["number cars"]
    nodal_transport_data.loc[nodal_transport_data["average fuel efficiency"] == 0.,"average fuel efficiency"] = transport_data["average fuel efficiency"].mean()


    #electric motors are more efficient, so alter transport demand

    #kWh/km from EPA https://www.fueleconomy.gov/feg/ for Tesla Model S
    plug_to_wheels_eta = 0.20
    battery_to_wheels_eta = plug_to_wheels_eta*0.9

    efficiency_gain = nodal_transport_data["average fuel efficiency"]/battery_to_wheels_eta


    #get heating demand for correction to demand time series
    temperature = xr.open_dataarray(snakemake.input.temp_air_total).T.to_pandas()

    #correction factors for vehicle heating
    dd_ICE = transport_degree_factor(temperature,
                                     options['transport_heating_deadband_lower'],
                                     options['transport_heating_deadband_upper'],
                                     options['ICE_lower_degree_factor'],
                                     options['ICE_upper_degree_factor'])

    dd_EV = transport_degree_factor(temperature,
                                    options['transport_heating_deadband_lower'],
                                    options['transport_heating_deadband_upper'],
                                    options['EV_lower_degree_factor'],
                                    options['EV_upper_degree_factor'])

    #divide out the heating/cooling demand from ICE totals
    ICE_correction = (transport_shape*(1+dd_ICE)).sum()/transport_shape.sum()

    transport = (transport_shape.multiply(nodal_energy_totals["total road"] + nodal_energy_totals["total rail"]
                                         - nodal_energy_totals["electricity rail"])*1e6*Nyears).divide(efficiency_gain*ICE_correction)

    #multiply back in the heating/cooling demand for EVs
    transport = transport.multiply(1+dd_EV)

    #Add transport demand factor depending on the year
    transport = transport * get_parameter(options["land_transport_demand"])

    ## derive plugged-in availability for PKW's (cars)


    traffic = pd.read_csv(os.path.join(snakemake.input.traffic_data, "Pkw__count"),
                          skiprows=2)["count"]

    avail_max = 0.95

    avail_mean = 0.8

    avail = avail_max - (avail_max - avail_mean)*(traffic - traffic.min())/(traffic.mean() - traffic.min())

    avail_profile = generate_periodic_profiles(dt_index=network.snapshots.tz_localize("UTC"),
                                               nodes=pop_layout.index,
                                               weekly_profile=avail.values).tz_localize(None)

    dsm_week = np.zeros((24*7,))

    dsm_week[(np.arange(0,7,1)*24+options['bev_dsm_restriction_time'])] = options['bev_dsm_restriction_value']

    dsm_profile = generate_periodic_profiles(dt_index=network.snapshots.tz_localize("UTC"),
                                             nodes=pop_layout.index,
                                             weekly_profile=dsm_week).tz_localize(None)


    ###############
    #CO2
    ###############

    #1e6 to convert Mt to tCO2
    co2_totals = 1e6*pd.read_csv(snakemake.input.co2_totals_name,index_col=0)



    return nodal_energy_totals, heat_demand, ashp_cop, gshp_cop, solar_thermal, transport, avail_profile, dsm_profile, co2_totals, nodal_transport_data



def prepare_costs(cost_file, USD_to_EUR, discount_rate, Nyears, lifetime):

    #set all asset costs and other parameters
    costs = pd.read_csv(cost_file,index_col=list(range(2))).sort_index()

    #correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"),"value"]*=1e3
    costs.loc[costs.unit.str.contains("USD"),"value"]*=USD_to_EUR

    #min_count=1 is important to generate NaNs which are then filled by fillna
    costs = costs.loc[:, "value"].unstack(level=1).groupby("technology").sum(min_count=1)
    costs = costs.fillna({"CO2 intensity" : 0,
                          "FOM" : 0,
                          "VOM" : 0,
                          "discount rate" : discount_rate,
                          "efficiency" : 1,
                          "fuel" : 0,
                          "investment" : 0,
                          "lifetime" : lifetime
    })

    costs["fixed"] = [(annuity(v["lifetime"],v["discount rate"])+v["FOM"]/100.)*v["investment"]*Nyears for i,v in costs.iterrows()]
    return costs

def add_generation(network):
    print("adding electricity generation")
    nodes = pop_layout.index

    conventionals = [("OCGT","gas")]

    for generator,carrier in [("OCGT","gas")]:
        network.add("Carrier",
                    carrier)

        network.madd("Bus",
                     ["EU " + carrier],
                     location="EU",
                     carrier=carrier)

        #use madd to get carrier inserted
        network.madd("Store",
                     ["EU " + carrier + " Store"],
                     bus=["EU " + carrier],
                     e_nom_extendable=True,
                     e_cyclic=True,
                     carrier=carrier,
                     capital_cost=0.) #could correct to e.g. 0.2 EUR/kWh * annuity and O&M

        network.add("Generator",
                    "EU " + carrier,
                    bus="EU " + carrier,
                    p_nom_extendable=True,
                    carrier=carrier,
                    capital_cost=0.,
                    marginal_cost=costs.at[carrier,'fuel'])


        network.madd("Link",
                     nodes + " " + generator,
                     bus0=["EU " + carrier]*len(nodes),
                     bus1=nodes,
                     bus2="co2 atmosphere",
                     marginal_cost=costs.at[generator,'efficiency']*costs.at[generator,'VOM'], #NB: VOM is per MWel
                     capital_cost=costs.at[generator,'efficiency']*costs.at[generator,'fixed'], #NB: fixed cost is per MWel
                     p_nom_extendable=True,
                     carrier=generator,
                     efficiency=costs.at[generator,'efficiency'],
                     efficiency2=costs.at[carrier,'CO2 intensity'],
                     lifetime=costs.at[generator,'lifetime'])

def add_wave(network, wave_cost_factor):
    wave_fn = "data/WindWaveWEC_GLTB.xlsx"

    locations = ["FirthForth","Hebrides"]

    #in kW
    capacity = pd.Series([750,1000,600],["Attenuator","F2HB","MultiPA"])

    #in EUR/MW
    costs = wave_cost_factor*pd.Series([2.5,2,1.5],["Attenuator","F2HB","MultiPA"])*1e6

    sheets = {}

    for l in locations:
        sheets[l] = pd.read_excel(wave_fn,
                                  index_col=0,skiprows=[0],parse_dates=True,
                                  sheet_name=l)

    to_drop = ["Vestas 3MW","Vestas 8MW"]
    wave = pd.concat([sheets[l].drop(to_drop,axis=1).divide(capacity,axis=1) for l in locations],
                     keys=locations,
                     axis=1)

    for wave_type in costs.index:
        n.add("Generator",
              "Hebrides "+wave_type,
              bus="GB4 0",
              p_nom_extendable=True,
              carrier="wave",
              capital_cost=(annuity(25,0.07)+0.03)*costs[wave_type],
              p_max_pu=wave["Hebrides",wave_type])



def insert_electricity_distribution_grid(network):
    print("Inserting electricity distribution grid with investment cost factor of",
          snakemake.config["sector"]['electricity_distribution_grid_cost_factor'])

    nodes = pop_layout.index

    network.madd("Bus",
                 nodes+ " low voltage",
                 location=nodes,
                 carrier="low voltage")

    network.madd("Link",
                 nodes + " electricity distribution grid",
                 bus0=nodes,
                 bus1=nodes + " low voltage",
                 p_nom_extendable=True,
                 p_min_pu=-1,
                 carrier="electricity distribution grid",
                 efficiency=1,
                 marginal_cost=0,
                 lifetime=costs.at['electricity distribution grid','lifetime'],
                 capital_cost=costs.at['electricity distribution grid','fixed']*snakemake.config["sector"]['electricity_distribution_grid_cost_factor'])


    #this catches regular electricity load and "industry electricity"
    loads = network.loads.index[network.loads.carrier.str.contains("electricity")]
    network.loads.loc[loads,"bus"] += " low voltage"

    bevs = network.links.index[network.links.carrier == "BEV charger"]
    network.links.loc[bevs,"bus0"] += " low voltage"

    v2gs = network.links.index[network.links.carrier == "V2G"]
    network.links.loc[v2gs,"bus1"] += " low voltage"

    hps = network.links.index[network.links.carrier.str.contains("heat pump")]
    network.links.loc[hps,"bus0"] += " low voltage"

    rh = network.links.index[network.links.carrier.str.contains("resistive heater")]
    network.links.loc[rh, "bus0"] += " low voltage"

    mchp = network.links.index[network.links.carrier.str.contains("micro gas")]
    network.links.loc[mchp, "bus1"] += " low voltage"

    #set existing solar to cost of utility cost rather the 50-50 rooftop-utility
    solar = network.generators.index[network.generators.carrier == "solar"]
    network.generators.loc[solar, "capital_cost"] = costs.at['solar-utility',
                                                             'fixed']
    if snakemake.wildcards.clusters[-1:] == "m":
        pop_solar = simplified_pop_layout.total.rename(index = lambda x: x + " solar")
    else:
        pop_solar = pop_layout.total.rename(index = lambda x: x + " solar")

    # add max solar rooftop potential assuming 0.1 kW/m2 and 10 m2/person,
    #i.e. 1 kW/person (population data is in thousands of people) so we get MW
    potential = 0.1*10*pop_solar

    network.madd("Generator",
                 solar,
                 suffix=" rooftop",
                 bus=network.generators.loc[solar, "bus"] + " low voltage",
                 carrier="solar rooftop",
                 p_nom_extendable=True,
                 p_nom_max=potential,
                 marginal_cost=network.generators.loc[solar, 'marginal_cost'],
                 capital_cost=costs.at['solar-rooftop', 'fixed'],
                 efficiency=network.generators.loc[solar, 'efficiency'],
                 p_max_pu=network.generators_t.p_max_pu[solar])


    network.add("Carrier","home battery")

    network.madd("Bus",
                 nodes + " home battery",
                 location=nodes,
                 carrier="home battery")

    network.madd("Store",
                 nodes + " home battery",
                 bus=nodes + " home battery",
                 e_cyclic=True,
                 e_nom_extendable=True,
                 carrier="home battery",
                 capital_cost=costs.at['battery storage','fixed'],
                 lifetime=costs.at['battery storage','lifetime'])

    network.madd("Link",
                 nodes + " home battery charger",
                 bus0=nodes + " low voltage",
                 bus1=nodes + " home battery",
                 carrier="home battery charger",
                 efficiency=costs.at['battery inverter','efficiency']**0.5,
                 capital_cost=costs.at['battery inverter','fixed'],
                 p_nom_extendable=True,
                 lifetime=costs.at['battery inverter','lifetime'])

    network.madd("Link",
                 nodes + " home battery discharger",
                 bus0=nodes + " home battery",
                 bus1=nodes + " low voltage",
                 carrier="home battery discharger",
                 efficiency=costs.at['battery inverter','efficiency']**0.5,
                 marginal_cost=options['marginal_cost_storage'],
                 p_nom_extendable=True,
                 lifetime=costs.at['battery inverter','lifetime'])


def insert_gas_distribution_costs(network):
    f_costs = options['gas_distribution_grid_cost_factor']
    print("Inserting gas distribution grid with investment cost\
          factor of", f_costs)

    # gas boilers
    gas_b = network.links[network.links.carrier.str.contains("gas boiler") &
                          (~network.links.carrier.str.contains("urban central"))].index
    network.links.loc[gas_b, "capital_cost"] += costs.loc['electricity distribution grid']["fixed"] * f_costs
    # micro CHPs
    mchp = network.links.index[network.links.carrier.str.contains("micro gas")]
    network.links.loc[mchp,  "capital_cost"] += costs.loc['electricity distribution grid']["fixed"] * f_costs

    # TODO: Research methane grid costs and implement!
    # biogas
    biogas = network.links.index[network.links.carrier.str.contains("digestible biomass to gas")]
    network.links.loc[biogas,  "capital_cost"] += costs.loc['electricity distribution grid']["fixed"] * f_costs

    biosng = network.links.index[network.links.carrier.str.contains("solid biomass to gas")]
    network.links.loc[biosng,  "capital_cost"] += costs.loc['electricity distribution grid']["fixed"] * f_costs

    # lowT steam methane
    mchp = network.links.index[network.links.carrier.str.contains("methane for lowT industry")]
    network.links.loc[mchp,  "capital_cost"] += costs.loc['electricity distribution grid']["fixed"] * f_costs

def add_electricity_grid_connection(network):

    carriers = ["onwind","solar"]

    gens = network.generators.index[network.generators.carrier.isin(carriers)]

    network.generators.loc[gens,"capital_cost"] += costs.at['electricity grid connection','fixed']

def add_storage(network):
    print("adding electricity storage")
    nodes = pop_layout.index

    network.add("Carrier","H2")


    network.madd("Bus",
                 nodes+ " H2",
                 location=nodes,
                 carrier="H2")

    network.madd("Link",
                 nodes + " H2 Electrolysis",
                 bus1=nodes + " H2",
                 bus0=nodes,
                 p_nom_extendable=True,
                 carrier="H2 Electrolysis",
                 efficiency=costs.at["electrolysis","efficiency"],
                 capital_cost=costs.at["electrolysis","fixed"],
                 lifetime=costs.at['electrolysis','lifetime'])

    network.madd("Link",
                 nodes + " H2 Fuel Cell",
                 bus0=nodes + " H2",
                 bus1=nodes,
                 p_nom_extendable=True,
                 carrier ="H2 Fuel Cell",
                 efficiency=costs.at["fuel cell","efficiency"],
                 capital_cost=costs.at["fuel cell","fixed"]*costs.at["fuel cell","efficiency"], #NB: fixed cost is per MWel
                 lifetime=costs.at['fuel cell','lifetime'])

    cavern_nodes = pd.DataFrame()

    if options['hydrogen_underground_storage']:
         h2_salt_cavern_potential = pd.read_csv(snakemake.input.h2_cavern,
                                                index_col=0,squeeze=True)
         h2_cavern_ct = h2_salt_cavern_potential[~h2_salt_cavern_potential.isna()]
         cavern_nodes = pop_layout[pop_layout.ct.isin(h2_cavern_ct.index)]

         h2_capital_cost = costs.at["hydrogen storage underground", "fixed"]

         # assumptions: weight storage potential in a country by population
         # TODO: fix with real geographic potentials
         #convert TWh to MWh with 1e6
         h2_pot = h2_cavern_ct.loc[cavern_nodes.ct]
         h2_pot.index = cavern_nodes.index
         h2_pot = h2_pot * cavern_nodes.fraction * 1e6

         network.madd("Store",
                      cavern_nodes.index + " H2 Store",
                      bus=cavern_nodes.index + " H2",
                      e_nom_extendable=True,
                      e_nom_max=h2_pot.values,
                      e_cyclic=True,
                      carrier="H2 Store",
                      capital_cost=h2_capital_cost)

    # hydrogen stored overground
    h2_capital_cost = costs.at["hydrogen storage tank", "fixed"]
    nodes_overground = nodes.symmetric_difference(cavern_nodes.index)

    network.madd("Store",
                 nodes_overground + " H2 Store",
                 bus=nodes_overground + " H2",
                 e_nom_extendable=True,
                 e_cyclic=True,
                 carrier="H2 Store",
                 capital_cost=h2_capital_cost)

    h2_links = pd.DataFrame(columns=["bus0","bus1","length"])
    prefix = "H2 pipeline "
    connector = " -> "
    attrs = ["bus0","bus1","length"]

    candidates = pd.concat([network.lines[attrs],network.links.loc[network.links.carrier == "DC",attrs]],
                           keys=["lines","links"])

    for candidate in candidates.index:
        buses = [candidates.at[candidate,"bus0"],candidates.at[candidate,"bus1"]]
        buses.sort()
        name = prefix + buses[0] + connector + buses[1]
        if name not in h2_links.index:
            h2_links.at[name,"bus0"] = buses[0]
            h2_links.at[name,"bus1"] = buses[1]
            h2_links.at[name,"length"] = candidates.at[candidate,"length"]

    #TODO Add efficiency losses
    network.madd("Link",
                 h2_links.index,
                 bus0=h2_links.bus0.values + " H2",
                 bus1=h2_links.bus1.values + " H2",
                 p_min_pu=-1,
                 p_nom_extendable=True,
                 length=h2_links.length.values,
                 capital_cost=costs.at['H2 pipeline','fixed']*h2_links.length.values,
                 carrier="H2 pipeline",
                 lifetime=costs.at['H2 pipeline','lifetime'])


    network.add("Carrier","battery")

    network.madd("Bus",
                 nodes + " battery",
                 location=nodes,
                 carrier="battery")

    network.madd("Store",
                 nodes + " battery",
                 bus=nodes + " battery",
                 e_cyclic=True,
                 e_nom_extendable=True,
                 carrier="battery",
                 capital_cost=costs.at['battery storage','fixed'],
                 lifetime=costs.at['battery storage','lifetime'])

    network.madd("Link",
                 nodes + " battery charger",
                 bus0=nodes,
                 bus1=nodes + " battery",
                 carrier="battery charger",
                 efficiency=costs.at['battery inverter','efficiency']**0.5,
                 capital_cost=costs.at['battery inverter','fixed'],
                 p_nom_extendable=True,
                 lifetime=costs.at['battery inverter','lifetime'])

    network.madd("Link",
                 nodes + " battery discharger",
                 bus0=nodes + " battery",
                 bus1=nodes,
                 carrier="battery discharger",
                 efficiency=costs.at['battery inverter','efficiency']**0.5,
                 marginal_cost=options['marginal_cost_storage'],
                 p_nom_extendable=True,
                 lifetime=costs.at['battery inverter','lifetime'])

    if options['methanation']:
        network.madd("Link",
                     nodes + " Sabatier",
                     bus0=nodes+" H2",
                     bus1=["EU gas"]*len(nodes),
                     bus2="co2 stored",
                     # bus3="co2 atmosphere",
                     p_nom_extendable=True,
                     carrier="Sabatier",
                     efficiency=costs.at["methanation","efficiency"],
                     efficiency2=-costs.at['gas', 'CO2 intensity']*costs.at["methanation","efficiency"],#+costs.at['methanation', 'CO2 stored'] * costs.at['methanation', 'capture rate'], # (1-0.24)*0.2376, #delta between input CO2 and the CO2 byproduct, assumed to go back to the store. should be correct for eta=0.91 ##costs.at["methanation","efficiency"]*costs.at['gas','CO2 intensity'],
                     # efficiency3=costs.at['methanation', 'CO2 stored'] * (1 - costs.at['methanation', 'capture rate']),
                     capital_cost=costs.at["methanation","fixed"],
                     lifetime=costs.at['methanation','lifetime'])

    #TODO: Check helmeth mass balance
    if options['helmeth']:
        network.madd("Link",
                     nodes + " helmeth",
                     bus0=nodes,
                     bus1=["EU gas"]*len(nodes),
                     bus2="co2 stored",
                     carrier="helmeth",
                     p_nom_extendable=True,
                     efficiency=costs.at["helmeth","efficiency"],
                     efficiency2=-costs.at["helmeth","efficiency"]*costs.at['gas','CO2 intensity'],
                     capital_cost=costs.at["helmeth","fixed"],
                     lifetime=costs.at['helmeth','lifetime'])


    if options['SMR']:
        network.madd("Link",
                     nodes + " SMR CC",
                     bus0=["EU gas"]*len(nodes),
                     bus1=nodes+" H2",
                     bus2="co2 atmosphere",
                     bus3="co2 stored",
                     p_nom_extendable=True,
                     carrier="SMR CC",
                     efficiency=costs.at["SMR CC","efficiency"],
                     efficiency2=costs.at['gas','CO2 intensity']*(1-options["cc_fraction"]),
                     efficiency3=costs.at['gas','CO2 intensity']*options["cc_fraction"],
                     capital_cost=costs.at["SMR CC","fixed"],
                     lifetime=costs.at['SMR CC','lifetime'])

        network.madd("Link",
                     nodes + " SMR",
                     bus0=["EU gas"]*len(nodes),
                     bus1=nodes+" H2",
                     bus2="co2 atmosphere",
                     p_nom_extendable=True,
                     carrier="SMR",
                     efficiency=costs.at["SMR","efficiency"],
                     efficiency2=costs.at['gas','CO2 intensity'],
                     capital_cost=costs.at["SMR","fixed"],
                     lifetime=costs.at['SMR','lifetime'])


def add_land_transport(network):

    print("adding land transport")

    fuel_cell_share = get_parameter(options["land_transport_fuel_cell_share"])
    electric_share = get_parameter(options["land_transport_electric_share"])
    ice_share = 1 - fuel_cell_share - electric_share

    print("shares of FCEV, EV and ICEV are",
          fuel_cell_share,
          electric_share,
          ice_share)

    if ice_share < 0:
        print("Error, more FCEV and EV share than 1.")
        sys.exit()

    nodes = pop_layout.index


    if electric_share > 0:

        network.add("Carrier","Li ion")

        network.madd("Bus",
                     nodes,
                     location=nodes,
                     suffix=" EV battery",
                     carrier="Li ion")

        network.madd("Load",
                     nodes,
                     suffix=" land transport EV",
                     bus=nodes + " EV battery",
                     carrier="land transport EV",
                     p_set=electric_share*(transport[nodes]+shift_df(transport[nodes],1)+shift_df(transport[nodes],2))/3.)

        p_nom = nodal_transport_data["number cars"]*0.011*electric_share  #3-phase charger with 11 kW * x% of time grid-connected

        network.madd("Link",
                     nodes,
                     suffix= " BEV charger",
                     bus0=nodes,
                     bus1=nodes + " EV battery",
                     p_nom=p_nom,
                     carrier="BEV charger",
                     p_max_pu=avail_profile[nodes],
                     efficiency=0.9, #[B]
                     #These were set non-zero to find LU infeasibility when availability = 0.25
                     #p_nom_extendable=True,
                     #p_nom_min=p_nom,
                     #capital_cost=1e6,  #i.e. so high it only gets built where necessary
        )

        if options["v2g"]:

            network.madd("Link",
                         nodes,
                         suffix=" V2G",
                         bus1=nodes,
                         bus0=nodes + " EV battery",
                         p_nom=p_nom,
                         carrier="V2G",
                         p_max_pu=avail_profile[nodes],
                         efficiency=0.9)  #[B]



        if options["bev_dsm"]:

            network.madd("Store",
                         nodes,
                         suffix=" battery storage",
                         bus=nodes + " EV battery",
                         carrier="battery storage",
                         e_cyclic=True,
                         e_nom=nodal_transport_data["number cars"]*0.05*options["bev_availability"]*electric_share, #50 kWh battery http://www.zeit.de/mobilitaet/2014-10/auto-fahrzeug-bestand
                         e_max_pu=1,
                         e_min_pu=dsm_profile[nodes])


    if fuel_cell_share > 0:

        network.madd("Load",
                     nodes,
                     suffix=" land transport fuel cell",
                     bus=nodes + " H2",
                     carrier="land transport fuel cell",
                     p_set=fuel_cell_share/options['transport_fuel_cell_efficiency']*transport[nodes])


    if ice_share > 0:

        if "EU oil" not in network.buses.index:
            network.madd("Bus",
                         ["EU oil"],
                         location="EU",
                         carrier="oil")

        network.madd("Load",
                     nodes,
                     suffix=" land transport oil",
                     bus="EU oil",
                     carrier="land transport oil",
                     p_set=ice_share/options['transport_internal_combustion_efficiency']*transport[nodes])

        # TODO: make hourly for each node? - the method now is inexact when running coarse resolution (e.g. > 1h)
        # co2 = ice_share/options['transport_internal_combustion_efficiency']*transport[nodes].sum().sum()/8760.*costs.at["oil",'CO2 intensity']
        network.madd("Load",
                     nodes,
                     suffix=" land transport oil emissions",
                     bus="co2 atmosphere",
                     carrier="land transport oil emissions",
                     p_set=-ice_share/options['transport_internal_combustion_efficiency']*transport[nodes]*costs.at["oil",'CO2 intensity']) #co2_4)

def add_heat(network):

    print("adding heat")

    sectors = ["residential", "services"]

    nodes = create_nodes_for_heat_sector()

    #NB: must add costs of central heating afterwards (EUR 400 / kWpeak, 50a, 1% FOM from Fraunhofer ISE)

    urban_fraction = options['central_fraction']*pop_layout["urban"]/(pop_layout[["urban","rural"]].sum(axis=1))

    # exogenously reduce space heat demand
    if options["reduce_space_heat_exogenously"]:
        dE = get_parameter(options["reduce_space_heat_exogenously_factor"])
        print("assumed space heat reduction of {} %".format(dE*100))
        for sector in sectors:
            heat_demand[sector + " space"] = (1-dE)*heat_demand[sector + " space"]

    heat_systems = ["residential rural", "services rural",
                    "residential urban decentral","services urban decentral",
                    "urban central"]
    for name in heat_systems:

        name_type = "central" if name == "urban central" else "decentral"

        network.add("Carrier",name + " heat")

        network.madd("Bus",
                     nodes[name] + " " + name + " heat",
                     location=nodes[name],
                     carrier=name + " heat")

        ## Add heat load

        for sector in sectors:
            if "rural" in name:
                factor = 1-urban_fraction[nodes[name]]
            elif "urban" in name:
                factor = urban_fraction[nodes[name]]
            else:
                factor = None
            if sector in name:
                heat_load = heat_demand[[sector + " water",sector + " space"]].groupby(level=1,axis=1).sum()[nodes[name]].multiply(factor)


        if name == "urban central":
            # TODO: seems to be an error in the grouping?
            heat_load = heat_demand.groupby(level=1,axis=1).sum()[nodes[name]].multiply(urban_fraction[nodes[name]]*(1+options['district_heating_loss']))
        network.madd("Load",
                     nodes[name],
                     suffix=" " + name + " heat",
                     bus=nodes[name] + " " + name + " heat",
                     carrier=name + " heat",
                     p_set=heat_load)


        ## Add heat pumps

        heat_pump_type = "air" if "urban" in name else "ground"

        costs_name = "{} {}-sourced heat pump".format(name_type,heat_pump_type)
        cop = {"air" : ashp_cop, "ground" : gshp_cop}
        efficiency = cop[heat_pump_type][nodes[name]] if options["time_dep_hp_cop"] else costs.at[costs_name,'efficiency']

        network.madd("Link",
                     nodes[name],
                     suffix=" {} {} heat pump".format(name,heat_pump_type),
                     bus0=nodes[name],
                     bus1=nodes[name] + " " + name + " heat",
                     carrier="{} {} heat pump".format(name,heat_pump_type),
                     efficiency=efficiency,
                     capital_cost=costs.at[costs_name,'efficiency']*costs.at[costs_name,'fixed'],
                     p_nom_extendable=True,
                     lifetime=costs.at[costs_name,'lifetime'])


        if options["tes"]:

            network.add("Carrier",name + " water tanks")

            network.madd("Bus",
                         nodes[name] + " " + name + " water tanks",
                         location=nodes[name],
                         carrier=name + " water tanks")

            network.madd("Link",
                         nodes[name] + " " + name + " water tanks charger",
                         bus0=nodes[name] + " " + name + " heat",
                         bus1=nodes[name] + " " + name + " water tanks",
                         efficiency=costs.at['water tank charger','efficiency'],
                         carrier=name + " water tanks charger",
                         p_nom_extendable=True)

            network.madd("Link",
                         nodes[name] + " " + name + " water tanks discharger",
                         bus0=nodes[name] + " " + name + " water tanks",
                         bus1=nodes[name] + " " + name + " heat",
                         carrier=name + " water tanks discharger",
                         efficiency=costs.at['water tank discharger','efficiency'],
                         p_nom_extendable=True)

            # [HP] 180 day time constant for centralised, 3 day for decentralised
            tes_time_constant_days = options["tes_tau"] if name_type == "decentral" else 180.

            network.madd("Store",
                         nodes[name] + " " + name + " water tanks",
                         bus=nodes[name] + " " + name + " water tanks",
                         e_cyclic=True,
                         e_nom_extendable=True,
                         carrier=name + " water tanks",
                         standing_loss=1-np.exp(-1/(24.*tes_time_constant_days)),
                         capital_cost=costs.at[name_type + ' water tank storage','fixed']/(1.17e-3*40), #conversion from EUR/m^3 to EUR/MWh for 40 K diff and 1.17 kWh/m^3/K
                         lifetime=costs.at[name_type + ' water tank storage','lifetime'])

        if options["boilers"]:

            network.madd("Link",
                         nodes[name] + " " + name + " resistive heater",
                         bus0=nodes[name],
                         bus1=nodes[name] + " " + name + " heat",
                         carrier=name + " resistive heater",
                         efficiency=costs.at[name_type + ' resistive heater','efficiency'],
                         capital_cost=costs.at[name_type + ' resistive heater','efficiency']*costs.at[name_type + ' resistive heater','fixed'],
                         p_nom_extendable=True,
                         lifetime=costs.at[name_type + ' resistive heater','lifetime'])

            network.madd("Link",
                         nodes[name] + " " + name + " gas boiler",
                         p_nom_extendable=True,
                         bus0=["EU gas"]*len(nodes[name]),
                         bus1=nodes[name] + " " + name + " heat",
                         bus2="co2 atmosphere",
                         carrier=name + " gas boiler",
                         efficiency=costs.at[name_type + ' gas boiler','efficiency'],
                         efficiency2=costs.at['gas','CO2 intensity'],
                         capital_cost=costs.at[name_type + ' gas boiler','efficiency']*costs.at[name_type + ' gas boiler','fixed'],
                         lifetime=costs.at[name_type + ' gas boiler','lifetime'])



        if options["solar_thermal"]:

            network.add("Carrier",name + " solar thermal")

            network.madd("Generator",
                         nodes[name],
                         suffix=" " + name + " solar thermal collector",
                         bus=nodes[name] + " " + name + " heat",
                         carrier=name + " solar thermal",
                         p_nom_extendable=True,
                         capital_cost=costs.at[name_type + ' solar thermal','fixed'],
                         p_max_pu=solar_thermal[nodes[name]],
                         lifetime=costs.at[name_type + ' solar thermal','lifetime'])


        if options["chp"]:

            if name == "urban central":
                #add gas CHP; biomass CHP is added in biomass section
                network.madd("Link",
                             nodes[name] + " urban central gas CHP",
                             bus0="EU gas",
                             bus1=nodes[name],
                             bus2=nodes[name] + " urban central heat",
                             bus3="co2 atmosphere",
                             carrier="urban central gas CHP",
                             p_nom_extendable=True,
                             capital_cost=costs.at['central gas CHP','fixed']*costs.at['central gas CHP','efficiency'],
                             marginal_cost=costs.at['central gas CHP','VOM'],
                             efficiency=costs.at['central gas CHP','efficiency'],
                             efficiency2=costs.at['central gas CHP','efficiency']/costs.at['central gas CHP','c_b'],
                             efficiency3=costs.at['gas','CO2 intensity'],
                             lifetime=costs.at['central gas CHP','lifetime'])

                network.madd("Link",
                             nodes[name] + " urban central gas CHP CC",
                             bus0="EU gas",
                             bus1=nodes[name],
                             bus2=nodes[name] + " urban central heat",
                             bus3="co2 atmosphere",
                             bus4="co2 stored",
                             carrier="urban central gas CHP CC",
                             p_nom_extendable=True,
                             capital_cost=costs.at['central gas CHP','fixed']*costs.at['central gas CHP','efficiency'] + costs.at['biomass CHP capture','fixed']*costs.at['gas','CO2 intensity'],
                             marginal_cost=costs.at['central gas CHP','VOM'],
                             efficiency=costs.at['central gas CHP','efficiency'] - costs.at['gas','CO2 intensity']*(costs.at['biomass CHP capture','electricity-input'] + costs.at['biomass CHP capture','compression-electricity-input']),
                             efficiency2=costs.at['central gas CHP','efficiency']/costs.at['central gas CHP','c_b'] + costs.at['gas','CO2 intensity']*(costs.at['biomass CHP capture','heat-output'] + costs.at['biomass CHP capture','compression-heat-output'] - costs.at['biomass CHP capture','heat-input']),
                             efficiency3=costs.at['gas','CO2 intensity']*(1-costs.at['biomass CHP capture','capture_rate']),
                             efficiency4=costs.at['gas','CO2 intensity']*costs.at['biomass CHP capture','capture_rate'],
                             lifetime=costs.at['central gas CHP','lifetime'])

            else:
                if options["micro_chp"]:
                    network.madd("Link",
                             nodes[name] + " " + name + " micro gas CHP",
                             p_nom_extendable=True,
                             bus0="EU gas",
                             bus1=nodes[name],
                             bus2=nodes[name] + " " + name + " heat",
                             bus3="co2 atmosphere",
                             carrier=name + " micro gas CHP",
                             efficiency=costs.at['micro CHP','efficiency'],
                             efficiency2=costs.at['micro CHP','efficiency-heat'],
                             efficiency3=costs.at['gas','CO2 intensity'],
                             capital_cost=costs.at['micro CHP','fixed'],
                             lifetime=costs.at['micro CHP','lifetime'])


    if options['retrofitting']['retro_endogen']:

        print("adding retrofitting endogenously")

        # resample heat demand temporal 'heat_demand_r' depending on in config
        # specified temporal resolution, to not overestimate retrofitting
        hours = list(filter(re.compile(r'^\d+h$', re.IGNORECASE).search, opts))
        if len(hours)==0:
            hours = [n.snapshots[1] - n.snapshots[0]]
        heat_demand_r =  heat_demand.resample(hours[0]).mean()

        # retrofitting data 'retro_data' with 'costs' [EUR/m^2] and heat
        # demand 'dE' [per unit of original heat demand] for each country and
        # different retrofitting strengths [additional insulation thickness in m]
        retro_data = pd.read_csv(snakemake.input.retro_cost_energy,
                                 index_col=[0, 1], skipinitialspace=True,
                                 header=[0, 1])
        # heated floor area [10^6 * m^2] per country
        floor_area = pd.read_csv(snakemake.input.floor_area, index_col=[0, 1])

        network.add("Carrier", "retrofitting")

        # share of space heat demand 'w_space' of total heat demand
        w_space = {}
        for sector in sectors:
            w_space[sector] = heat_demand_r[sector + " space"] / \
                (heat_demand_r[sector + " space"] + heat_demand_r[sector + " water"])
        w_space["tot"] = ((heat_demand_r["services space"] +
                           heat_demand_r["residential space"]) /
                           heat_demand_r.groupby(level=[1], axis=1).sum())


        for name in network.loads[network.loads.carrier.isin([x + " heat" for x in heat_systems])].index:

            node = network.buses.loc[name, "location"]
            ct = pop_layout.loc[node, "ct"]

            # weighting 'f' depending on the size of the population at the node
            f = urban_fraction[node] if "urban" in name else (1-urban_fraction[node])
            if f == 0:
                continue
            # get sector name ("residential"/"services"/or both "tot" for urban central)
            sec = [x if x in name else "tot" for x in sectors][0]

            # get floor aread at node and region (urban/rural) in m^2
            floor_area_node = ((pop_layout.loc[node].fraction
                                  * floor_area.loc[ct, "value"] * 10**6).loc[sec] * f)
            # total heat demand at node [MWh]
            demand = (network.loads_t.p_set[name].resample(hours[0])
                      .mean())

            # space heat demand at node [MWh]
            space_heat_demand = demand * w_space[sec][node]
            # normed time profile of space heat demand 'space_pu' (values between 0-1),
            # p_max_pu/p_min_pu of retrofitting generators
            space_pu = (space_heat_demand / space_heat_demand.max()).to_frame(name=node)

            # minimum heat demand 'dE' after retrofitting in units of original heat demand (values between 0-1)
            dE = retro_data.loc[(ct, sec), ("dE")]
            # get addtional energy savings 'dE_diff' between the different retrofitting strengths/generators at one node
            dE_diff = abs(dE.diff()).fillna(1-dE.iloc[0])
            # convert costs Euro/m^2 -> Euro/MWh
            capital_cost =  retro_data.loc[(ct, sec), ("cost")] * floor_area_node / \
                            ((1 - dE) * space_heat_demand.max())
            # number of possible retrofitting measures 'strengths' (set in list at config.yaml 'l_strength')
            # given in additional insulation thickness [m]
            # for each measure, a retrofitting generator is added at the node
            strengths = retro_data.columns.levels[1]

            # check that ambitious retrofitting has higher costs per MWh than moderate retrofitting
            if (capital_cost.diff() < 0).sum():
                print(
                    "warning, costs are not linear for ", ct, " ", sec)
                s = capital_cost[(capital_cost.diff() < 0)].index
                strengths = strengths.drop(s)

            # reindex normed time profile of space heat demand back to hourly resolution
            space_pu = (space_pu.reindex(index=heat_demand.index)
                          .fillna(method="ffill"))

            # add for each retrofitting strength a generator with heat generation profile following the profile of the heat demand
            for strength in strengths:
                network.madd('Generator',
                            [node],
                            suffix=' retrofitting ' + strength + " " + name[6::],
                            bus=name,
                            carrier="retrofitting",
                            p_nom_extendable=True,
                            p_nom_max=dE_diff[strength] * space_heat_demand.max(), # maximum energy savings for this renovation strength
                            p_max_pu=space_pu,
                            p_min_pu=space_pu,
                            country=ct,
                            capital_cost=capital_cost[strength] * options['retrofitting']['cost_factor'])


def create_nodes_for_heat_sector():
    sectors = ["residential", "services"]
    # stores the different groups of nodes
    nodes = {}
    # rural are areas with low heating density and individual heating
    # urban are areas with high heating density
    # urban can be split into district heating (central) and individual heating (decentral)
    for sector in sectors:
        nodes[sector + " rural"] = pop_layout.index

        if options["central"]:
            urban_decentral_ct = pd.Index(["ES", "GR", "PT", "IT", "BG"])
            nodes[sector + " urban decentral"] = pop_layout.index[pop_layout.ct.isin(urban_decentral_ct)]
        else:
            nodes[sector + " urban decentral"] = pop_layout.index
    # for central nodes, residential and services are aggregated
    nodes["urban central"] = pop_layout.index.symmetric_difference(nodes["residential urban decentral"])
    return nodes


def add_biomass(network):

    print("adding biomass")

    nodes = pop_layout.index

    #biomass distributed at country level - i.e. transport within country allowed
    cts = pop_layout.ct.value_counts().index

    biomass_potentials = pd.read_csv(snakemake.input.biomass_potentials,
                                     index_col=0)

    # potential per node distributed within country by population
    biomass_pot_node = (biomass_potentials.loc[pop_layout.ct]
                        .set_index(pop_layout.index)
                        .mul(pop_layout.fraction, axis="index"))

    print(biomass_pot_node)
    network.add("Carrier","digestible biomass")

    network.madd("Bus",
                 nodes + " digestible biomass",
                 location=nodes,
                 carrier="digestible biomass")

    network.madd("Store",
                 nodes + " digestible biomass",
                 bus=nodes + " digestible biomass",
                 carrier="digestible biomass",
                 e_cyclic=True)

    network.add("Carrier","solid biomass")

    network.madd("Bus",
                 nodes + " solid biomass",
                 location=nodes,
                 carrier="solid biomass")

    network.madd("Store",
                 nodes + " solid biomass",
                 bus=nodes + " solid biomass",
                 carrier="solid biomass",
                 e_cyclic=True)

    digestible_biomass_types = ["manureslurry", "municipal biowaste",
                                "sewage sludge", "straw"] #"Food industry residues",
                                # "Landscape management", "Straw"]

    biomass_potential = {}
    biomass_cost = {}

    for name in digestible_biomass_types:

        biomass_potential[name] = biomass_pot_node[name].values
        network.add("Carrier", name + " digestible biomass")

        network.madd("Bus",
                 nodes + " " + name + " digestible biomass",
                 location=nodes,
                 carrier=name + " digestible biomass")

        #TODO: add individual costs
        # network.madd("Store",
        #              nodes + " " + name + " digestible biomass",
        #              bus=nodes + " " + name + " digestible biomass",
        #              carrier=name + " digestible biomass",
        #              e_nom=biomass_pot_node[name].values,
        #              # marginal_cost=costs.at['digestible biomass', 'fuel'], #Add individual costs
        #              e_initial=biomass_pot_node[name].values)

        network.madd("Store",
                    nodes + " " + name + " digestible biomass",
                    bus=nodes + " " + name + " digestible biomass",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    carrier=name + " digestible biomass")

        network.madd("Generator",
                     nodes + " " + name + " digestible biomass",
                     bus=nodes + " " + name + " digestible biomass",
                     carrier=name + " digestible biomass",
                     p_nom_extendable=True,
                     p_nom_max=biomass_potential[name]/8760,
                     marginal_cost=costs.at['digestible biomass', 'fuel'])


        #NB: Assuming that the input substrates are given in the biogas potential and biogas cost!

        network.madd("Link",
                     nodes + " " + name + " digestible biomass",
                     bus0=nodes + " " + name + " digestible biomass",
                     bus1=nodes + " digestible biomass",
                     bus2="co2 atmosphere",
                     carrier="digestible biomass",
                     efficiency=1.,
                     efficiency2=-costs.at['gas', 'CO2 intensity'] - costs.at["Anaerobic digestion", "CO2 stored"],  #Adding the CO2 in the biogas mix
                     p_nom_extendable=True)

    #TODO: gas grid cost added for biogas processes in insert_gas_distribution_costs, but demands refining! Also add CO2 transport cost!
    # network.madd("Link",
    #              nodes + " digestible biomass",
    #              bus0=nodes + " digestible biomass",
    #              bus1="EU gas",
    #              bus3="co2 atmosphere",
    #              carrier="digestible biomass to gas",
    #              capital_cost=costs.at["Anaerobic digestion", "fixed"] + costs.at["biogas upgrading", "fixed"],  # Change to DEA values
    #              marginal_cost=costs.at["biogas upgrading", "VOM"],  # Add DEA values for biogas process
    #              efficiency=1,
    #              efficiency3=costs.at["Anaerobic digestion", "CO2 stored"],
    #              p_nom_extendable=True)

    network.madd("Link",
                 nodes + " digestible biomass CC",
                 bus0=nodes + " digestible biomass",
                 bus1="EU gas",
                 bus2="co2 stored",
                 bus3="co2 atmosphere",
                 carrier="digestible biomass to gas",
                 capital_cost=costs.at["Anaerobic digestion", "fixed"] + costs.at["biogas upgrading", "fixed"], #Assuming that the CO2 from upgrading is pure, such as in amine scrubbing. I.e., with and without CC is equivalent.
                 marginal_cost=costs.at["biogas upgrading", "VOM"],  # Add DEA values for biogas process
                 efficiency=1,
                 efficiency2=costs.at["Anaerobic digestion", "CO2 stored"] * costs.at['Anaerobic digestion','capture rate'],
                 efficiency3=costs.at["Anaerobic digestion", "CO2 stored"] * (1-costs.at['Anaerobic digestion','capture rate']),
                 p_nom_extendable=True)

    network.madd("Link",
                 nodes + " digestible biomass to hydrogen CC",
                 bus0=nodes + " digestible biomass",
                 bus1=nodes + " H2",
                 bus2="co2 stored",
                 bus3="co2 atmosphere",
                 carrier="digestible biomass to hydrogen",
                 capital_cost=costs.at['BtL', 'fixed'] + costs.at["biogas upgrading", "fixed"],
                 marginal_cost=costs.at["biogas upgrading", "VOM"],
                 efficiency=0.45,
                 efficiency2=(costs.at['gas', 'CO2 intensity'] + costs.at["Anaerobic digestion", "CO2 stored"]) * costs.at['Anaerobic digestion','capture rate'],
                 efficiency3=(costs.at['gas', 'CO2 intensity'] + costs.at["Anaerobic digestion", "CO2 stored"]) * (1-costs.at['Anaerobic digestion','capture rate']),
                 p_nom_extendable=True)


    solid_biomass_types = ["poplar", "forest residues", "industry wood residues"] #"scrap wood"


    for name in solid_biomass_types:
        network.add("Carrier", name + " solid biomass")

        network.madd("Bus",
                     nodes + " " + name + " solid biomass",
                     location=nodes,
                     carrier=name + " solid biomass")

        biomass_potential[name] = biomass_pot_node[name].values
        #TODO: change to individual biomass type cost, per node
        biomass_cost[name] = costs.at['solid biomass', 'fuel']

        network.madd("Store",
                    nodes + " " + name + " solid biomass",
                    bus=nodes + " " + name + " solid biomass",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    carrier=name + " solid biomass")

        network.madd("Generator",
                     nodes + " " + name + " solid biomass",
                     bus=nodes + " " + name + " solid biomass",
                     carrier=name + " solid biomass",
                     p_nom_extendable=True,
                     p_nom_max=biomass_potential[name]/8760,
                     marginal_cost=biomass_cost[name])

        network.madd("Link",
                     nodes + " " + name + " solid biomass",
                     bus0=nodes + " " + name + " solid biomass",
                     bus1=nodes + " solid biomass",
                     bus2="co2 atmosphere",
                     carrier="solid biomass",
                     efficiency=1.,
                     efficiency2=-costs.at['solid biomass', 'CO2 intensity'],
                     p_nom_extendable=True)


    # Add solid biomass import
    # if options["biomass_import"]:
    for o in opts:
        if o[o.find("B") + 4:] == "Im":
            print("Adding biomass import")

            network.add("Carrier","solid biomass import")
            import_potential = {}
            import_cost = {}
            import_name = {}
            superfluous = {}
            tot_EU_biomass = biomass_pot_node.values.sum() - biomass_pot_node["not included"].values.sum()
            print('Total EU biomass: ', tot_EU_biomass*3.6/1e9)
            step_size = 10 #EJ
            biomass_import_limit_low_level = 20e9 #EJ

            for num in range(1, 4):
                import_name[num] = "import" + str(num)
                if num == 1:
                    import_potential[num] = max(biomass_import_limit_low_level / 3.6 - tot_EU_biomass,0) # substract EU biomass from 20 EJ. If EU biomass > 20, return 0
                    import_cost[num] = 15 * 3.6  # EUR/MWh
                    superfluous = min(biomass_import_limit_low_level / 3.6 - tot_EU_biomass, 0) #If EU biomass > 20, reduce the following group(s)
                else:
                    import_potential[num] = max(step_size*1e9 / 3.6 + superfluous,0)  # EJ --> MWh
                    import_cost[num] = (15 + step_size*0.25 * (num - 1)) * 3.6  # EUR/MWh
                    superfluous += min(-superfluous, step_size*1e9 / 3.6)

                network.madd("Bus",
                            [import_name[num] + " solid biomass"],
                            location="EU",
                            carrier="solid biomass import")

                network.madd("Store",
                             [import_name[num] + " solid biomass"],
                             bus=import_name[num] + " solid biomass",
                             e_nom_extendable=True,
                             e_cyclic=True,
                             carrier="solid biomass import",
                             )

                network.madd("Generator",
                            [import_name[num] + " solid biomass"],
                            bus=import_name[num] + " solid biomass",
                            carrier="solid biomass import",
                            p_nom_extendable=True,
                            p_nom_max=import_potential[num]/8760,
                            marginal_cost=import_cost[num])

                network.madd("Link",
                            nodes + " " + import_name[num] + " solid biomass",
                            bus0=import_name[num] + " solid biomass",
                            bus1=nodes + " solid biomass",
                            bus2="co2 atmosphere",
                            carrier="solid biomass",
                            efficiency=1.,
                            efficiency2=-costs.at['solid biomass', 'CO2 intensity'], #Here, if the store at bus1 is cyclic
                            p_nom_extendable=True)


    #TODO: Add marginal costs
    network.madd("Link",
                 nodes + " solid biomass to gas",
                 bus0=nodes + " solid biomass",
                 bus1="EU gas",
                 bus3="co2 atmosphere",
                 carrier="BioSNG",
                 lifetime=costs.at['BioSNG', 'lifetime'],
                 efficiency=costs.at['BioSNG', 'efficiency'],
                 efficiency3=costs.at['BioSNG', 'CO2 stored'],
                 p_nom_extendable=True,
                 capital_cost=costs.at['BioSNG', 'fixed'],
                 marginal_cost=0., #costs.loc["BioSNG", "VOM"]
                 )

    network.madd("Link",
                 nodes + " solid biomass to gas CC",
                 bus0=nodes + " solid biomass",
                 bus1="EU gas",
                 bus2="co2 stored",
                 bus3="co2 atmosphere",
                 carrier="BioSNG",
                 lifetime=costs.at['BioSNG', 'lifetime'],
                 efficiency=costs.at['BioSNG', 'efficiency'],
                 efficiency2=costs.at['BioSNG', 'CO2 stored'] * costs.at['BioSNG', 'capture rate'],
                 efficiency3=costs.at['BioSNG', 'CO2 stored'] * (1 - costs.at['BioSNG', 'capture rate']),
                 p_nom_extendable=True,
                 capital_cost=costs.at['BioSNG', 'fixed'] + costs.at['biomass CHP capture','fixed']*costs.at["BioSNG", "CO2 stored"],
                 marginal_cost=0., #costs.loc["BioSNG", "VOM"]
                 )


    network.madd("Link",
                 nodes + " solid biomass to hydrogen",
                 bus0=nodes + " solid biomass",
                 bus1=nodes + " H2",
                 bus2="co2 stored",
                 bus3="co2 atmosphere",
                 carrier="solid biomass to hydrogen",
                 #lifetime=costs.at['BioSNG', 'lifetime'],
                 efficiency=0.5,#costs.at['BioSNG', 'efficiency'],
                 efficiency2=costs.at['solid biomass', 'CO2 intensity'] * costs.at['BtL', 'capture rate'],
                 efficiency3=costs.at['solid biomass', 'CO2 intensity'] * (1 - costs.at['BtL', 'capture rate']),
                 p_nom_extendable=True,
                 capital_cost=costs.at['BtL', 'fixed'], #CO2 separation included
                 marginal_cost=0., #costs.loc["BioSNG", "VOM"]
                 )

    #
    # network.madd("Bus",
    #              nodes + " liquid biofuel",
    #              location=nodes,
    #              carrier="biomass to liquid")
    #
    # network.madd("Store",
    #              nodes + " liquid biofuel",
    #              bus=nodes + " liquid biofuel",
    #              e_nom_extendable=True,
    #              e_cyclic=True,
    #              carrier="biomass to liquid")

    network.madd("Link",
                 nodes + " biomass to liquid",
                 bus0=nodes + " solid biomass",
                 bus1="EU oil",#nodes + " liquid biofuel",#"EU oil",
                 bus3="co2 atmosphere",
                 carrier="biomass to liquid",
                 lifetime=costs.at['BtL', 'lifetime'],
                 efficiency=costs.at['BtL', 'efficiency'],
                 efficiency3=costs.at['BtL', 'CO2 stored'],
                 p_nom_extendable=True,
                 capital_cost=costs.at['BtL', 'fixed'],
                 marginal_cost=0.,  # costs.loc["BtL", "VOM"]
                 )

    network.madd("Link",
                 nodes + " biomass to liquid CC",
                 bus0=nodes + " solid biomass",
                 bus1="EU oil",#nodes + " liquid biofuel",#"EU oil",
                 bus2="co2 stored",
                 bus3="co2 atmosphere",
                 carrier="biomass to liquid",
                 lifetime=costs.at['BtL', 'lifetime'],
                 efficiency=costs.at['BtL', 'efficiency'],
                 efficiency2=costs.at['BtL', 'CO2 stored'] * costs.at['BtL', 'capture rate'],
                 efficiency3=costs.at['BtL', 'CO2 stored'] * (1 - costs.at['BtL', 'capture rate']),
                 p_nom_extendable=True,
                 capital_cost=costs.at['BtL', 'fixed'] + costs.at['biomass CHP capture','fixed']*costs.at["BtL", "CO2 stored"],
                 marginal_cost=0.,  # costs.loc["BtL", "VOM"]
                 )
    #
    # network.madd("Link",
    #              nodes + " liquid biofuel",
    #              bus0=nodes + " liquid biofuel",
    #              bus1="EU oil",
    #              carrier="biomass to liquid",
    #              p_nom_extendable=True,
    #              )

    #TODO: Add real data for bioelectricity without CHP!
    network.madd("Link",
                 nodes + " solid biomass to electricity",
                 bus0=nodes + " solid biomass",
                 bus1=nodes,
                 bus3="co2 atmosphere",
                 carrier="solid biomass to electricity",
                 p_nom_extendable=True,
                 capital_cost=0.7 * costs.at['central solid biomass CHP', 'fixed'] * costs.at[
                     'central solid biomass CHP', 'efficiency'],
                 marginal_cost=costs.at['central solid biomass CHP', 'VOM'],
                 efficiency=0.4,
                 efficiency3=costs.at['solid biomass', 'CO2 intensity'],
                 lifetime=costs.at['central solid biomass CHP', 'lifetime'])

    network.madd("Link",
                 nodes + " solid biomass to electricity CC",
                 bus0=nodes + " solid biomass",
                 bus1=nodes,
                 bus3="co2 atmosphere",
                 bus4="co2 stored",
                 carrier="solid biomass to electricity CC",
                 p_nom_extendable=True,
                 capital_cost=0.7 * costs.at['central solid biomass CHP', 'fixed'] * costs.at[
                     'central solid biomass CHP', 'efficiency']
                              + costs.at['biomass CHP capture','fixed']*costs.at['solid biomass', 'CO2 intensity'],
                 marginal_cost=costs.at['central solid biomass CHP', 'VOM'],
                 efficiency=costs.at['central solid biomass CHP', 'efficiency'],
                 efficiency3=costs.at['solid biomass', 'CO2 intensity'] * (1 - options["cc_fraction"]),
                 efficiency4=costs.at['solid biomass', 'CO2 intensity'] * options["cc_fraction"],
                 # p_nom_ratio=costs.at['central solid biomass CHP', 'p_nom_ratio'],
                 lifetime=costs.at['central solid biomass CHP', 'lifetime'])

    #AC buses with district heating
    urban_central = network.buses.index[network.buses.carrier == "urban central heat"]
    if not urban_central.empty and options["chp"]:
        urban_central = urban_central.str[:-len(" urban central heat")]

        network.madd("Link",
                     urban_central + " urban central solid biomass CHP",
                     bus0=urban_central + " solid biomass",
                     bus1=urban_central,
                     bus2=urban_central + " urban central heat",
                     bus3="co2 atmosphere",
                     carrier="urban central solid biomass CHP",
                     p_nom_extendable=True,
                     capital_cost=costs.at['central solid biomass CHP','fixed']*costs.at['central solid biomass CHP','efficiency'],
                     marginal_cost=costs.at['central solid biomass CHP','VOM'],
                     efficiency=costs.at['central solid biomass CHP','efficiency'],
                     efficiency2=costs.at['central solid biomass CHP','efficiency-heat'],
                     efficiency3=costs.at['solid biomass', 'CO2 intensity'],
                     lifetime=costs.at['central solid biomass CHP','lifetime'])


        network.madd("Link",
                     urban_central + " urban central solid biomass CHP CC",
                     bus0=urban_central + " solid biomass",
                     bus1=urban_central,
                     bus2=urban_central + " urban central heat",
                     bus3="co2 atmosphere",
                     bus4="co2 stored",
                     carrier="urban central solid biomass CHP CC",
                     p_nom_extendable=True,
                     capital_cost=costs.at['central solid biomass CHP','fixed']*costs.at['central solid biomass CHP','efficiency']
                              + costs.at['biomass CHP capture','fixed']*costs.at['solid biomass', 'CO2 intensity'],
                     marginal_cost=costs.at['central solid biomass CHP','VOM'],
                     efficiency=costs.at['central solid biomass CHP','efficiency'],
                     efficiency2=costs.at['central solid biomass CHP', 'efficiency-heat'] + costs.at[
                         'solid biomass', 'CO2 intensity'] * (costs.at['biomass CHP capture', 'heat-output'] + costs.at[
                         'biomass CHP capture', 'compression-heat-output'] - costs.at[
                                                                  'biomass CHP capture', 'heat-output']),
                     efficiency3=costs.at['solid biomass','CO2 intensity']*(1-options["cc_fraction"]),
                     efficiency4=costs.at['solid biomass','CO2 intensity']*options["cc_fraction"],
                     c_b=costs.at['central solid biomass CHP','c_b'],
                     c_v=costs.at['central solid biomass CHP','c_v'],
                     p_nom_ratio=costs.at['central solid biomass CHP','p_nom_ratio'],
                     lifetime=costs.at['central solid biomass CHP','lifetime'])



def add_industry(network):

    print("adding industrial demand")

    nodes = pop_layout.index

    #1e6 to convert TWh to MWh
    industrial_demand = 1e6*pd.read_csv(snakemake.input.industrial_demand,
                                        index_col=0)

    network.madd("Bus",
                 nodes + " lowT process steam",
                 carrier="lowT process steam")


    network.madd("Load",
                 nodes,
                 suffix=" lowT process steam",
                 bus=nodes + " lowT process steam",
                 carrier="lowT process steam",
                 p_set=industrial_demand.loc[nodes,"solid biomass"]/8760.)

    for o in opts:
        if "B" in o:
            network.madd("Link",
                        nodes,
                        suffix=" solid biomass for lowT industry",
                        bus0=nodes + " solid biomass",
                        bus1=nodes + " lowT process steam",
                        bus2="co2 atmosphere",
                        carrier="lowT process steam solid biomass",
                        p_nom_extendable=True,
                        efficiency=costs.at['solid biomass to steam','efficiency'],
                        efficiency2=costs.at['solid biomass','CO2 intensity'],
                        capital_cost=costs.at['solid biomass to steam','fixed'])
            # TODO: cement capture rather low-cost compared to other processes (high conc. CO2) --> adapt!

            network.madd("Link",
                         nodes,
                         suffix=" solid biomass for lowT industry CC",
                         bus0=nodes + " solid biomass",
                         bus1=nodes + " lowT process steam",
                         bus2="co2 atmosphere",
                         bus3="co2 stored",
                         carrier="lowT process steam solid biomass CC",
                         p_nom_extendable=True,
                         efficiency=0.9 * costs.at['solid biomass to steam', 'efficiency'],
                         capital_cost=costs.at['solid biomass to steam', 'fixed'] + costs.at[
                             "cement capture", "fixed"] * costs.at['solid biomass', 'CO2 intensity'],
                         efficiency2=costs.at['solid biomass', 'CO2 intensity'] * (
                                     1 - costs.at["cement capture", "capture_rate"]),
                         efficiency3=costs.at['solid biomass', 'CO2 intensity'] * costs.at[
                             "cement capture", "capture_rate"],
                         lifetime=costs.at['cement capture', 'lifetime'])


    network.madd("Link",
                 nodes,
                 suffix=" methane for lowT industry",
                 bus0="EU gas",
                 bus1=nodes + " lowT process steam",
                 bus2="co2 atmosphere",
                 carrier="lowT process steam methane",
                 p_nom_extendable=True,
                 efficiency=costs.at['gas to steam','efficiency'],
                 capital_cost=costs.at['gas to steam','fixed'],
                 efficiency2=costs.at['gas','CO2 intensity'])

    network.madd("Link",
                 nodes,
                 suffix=" methane for lowT industry CC",
                 bus0="EU gas",
                 bus1=nodes + " lowT process steam",
                 bus2="co2 atmosphere",
                 bus3="co2 stored",
                 carrier="lowT process steam methane CC",
                 p_nom_extendable=True,
                 efficiency=0.9*costs.at['gas to steam','efficiency'],
                 capital_cost=costs.at['gas to steam','fixed']+costs.at["cement capture","fixed"]*costs.at['gas','CO2 intensity'],
                 efficiency2=costs.at['gas','CO2 intensity']*(1-costs.at["cement capture","capture_rate"]),
                 efficiency3=costs.at['gas','CO2 intensity']*costs.at["cement capture","capture_rate"],
                 lifetime=costs.at['cement capture','lifetime'])

    #TODO: adapt this to real values for H2 to steam
    network.madd("Link",
                 nodes,
                 suffix=" H2 for lowT industry",
                 bus0=nodes + " H2",
                 bus1=nodes + " lowT process steam",
                 carrier="lowT process steam H2",
                 p_nom_extendable=True,
                 efficiency=costs.at['gas to steam','efficiency'],
                 capital_cost=costs.at['gas to steam','fixed'])


    network.madd("Bus",
                 ["gas for industry"],
                 location="EU",
                 carrier="gas for industry")

    network.madd("Load",
                 ["gas for industry"],
                 bus="gas for industry",
                 carrier="gas for industry",
                 p_set=industrial_demand.loc[nodes,"methane"].sum()/8760.)

    network.madd("Link",
                 ["gas for industry"],
                 bus0="EU gas",
                 bus1="gas for industry",
                 bus2="co2 atmosphere",
                 carrier="gas for industry",
                 p_nom_extendable=True,
                 efficiency=1.,
                 efficiency2=costs.at['gas','CO2 intensity'])

    network.madd("Link",
                 ["gas for industry CC"],
                 bus0="EU gas",
                 bus1="gas for industry",
                 bus2="co2 atmosphere",
                 bus3="co2 stored",
                 carrier="gas for industry CC",
                 p_nom_extendable=True,
                 capital_cost=costs.at["cement capture","fixed"]*costs.at['gas','CO2 intensity'],
                 efficiency=0.9,
                 efficiency2=costs.at['gas','CO2 intensity']*(1-costs.at["cement capture","capture_rate"]),
                 efficiency3=costs.at['gas','CO2 intensity']*costs.at["cement capture","capture_rate"],
                 lifetime=costs.at['cement capture','lifetime'])


    network.madd("Load",
                 nodes,
                 suffix=" H2 for industry",
                 bus=nodes + " H2",
                 carrier="H2 for industry",
                 p_set=industrial_demand.loc[nodes,"hydrogen"]/8760.)


    network.madd("Load",
                 ["oil for shipping"],#nodes,
                 # suffix=" oil for shipping",
                 bus="EU oil",
                 carrier="oil for shipping",
                 p_set = (1-options['shipping_H2_share'])*get_parameter(options["shipping_demand"])*nodal_energy_totals.loc[nodes,["total international navigation","total domestic navigation"]].sum(axis=1).sum()*1e6/8760.)


    network.madd("Load",
                 nodes,
                 suffix=" H2 for shipping",
                 bus=nodes + " H2",
                 carrier="H2 for shipping",
                 p_set = options['shipping_H2_share']*get_parameter(options["shipping_demand"])*nodal_energy_totals.loc[nodes,["total international navigation","total domestic navigation"]].sum(axis=1)*1e6*options['shipping_average_efficiency']/costs.at["fuel cell","efficiency"]/8760.)


    if "EU oil" not in network.buses.index:
        network.madd("Bus",
                     ["EU oil"],
                     location="EU",
                     carrier="oil")

    #use madd to get carrier inserted
    if "EU oil Store" not in network.stores.index:
        network.madd("Store",
                     ["EU oil Store"],
                     bus="EU oil",
                     e_nom_extendable=True,
                     e_cyclic=True,
                     carrier="oil",
                     capital_cost=0.) #could correct to e.g. 0.001 EUR/kWh * annuity and O&M

    if "EU oil" not in network.generators.index:
        network.add("Generator",
                    "EU oil",
                    bus="EU oil",
                    p_nom_extendable=True,
                    carrier="oil",
                    capital_cost=0.,
                    marginal_cost=costs.at["oil",'fuel'])

    if options["oil_boilers"]:

        nodes_heat = create_nodes_for_heat_sector()

        for name in ["residential rural", "services rural", "residential urban decentral", "services urban decentral"]:
            network.madd("Link",
                         nodes_heat[name] + " " + name + " oil boiler",
                         p_nom_extendable=True,
                         bus0="EU oil",
                         bus1=nodes_heat[name] + " " + name + " heat",
                         bus2="co2 atmosphere",
                         carrier=name + " oil boiler",
                         efficiency=costs.at['decentral oil boiler', 'efficiency'],
                         efficiency2=costs.at['oil', 'CO2 intensity'],
                         capital_cost=costs.at['decentral oil boiler', 'efficiency'] * costs.at[
                                                'decentral oil boiler', 'fixed'],
                         lifetime=costs.at['decentral oil boiler','lifetime'])

    network.madd("Link",
                 nodes + " Fischer-Tropsch",
                 bus0=nodes + " H2",
                 bus1="EU oil",
                 bus2="co2 stored",
                 bus3="co2 atmosphere",
                 carrier="electrofuel",
                 efficiency=costs.at["Fischer-Tropsch",'efficiency'],
                 capital_cost=costs.at["Fischer-Tropsch",'fixed'],
                 efficiency2=-(1+(1-costs.at["Fischer-Tropsch",'capture rate']))*costs.at['oil', 'CO2 intensity']*costs.at["Fischer-Tropsch",'efficiency'],
                 efficiency3=(1-costs.at["Fischer-Tropsch",'capture rate'])*costs.at['oil', 'CO2 intensity']*costs.at["Fischer-Tropsch",'efficiency'],
                 p_nom_extendable=True,
                 lifetime=costs.at['Fischer-Tropsch','lifetime'])

    network.madd("Load",
                 ["naphtha for industry"],
                 bus="EU oil",
                 carrier="naphtha for industry",
                 p_set = industrial_demand.loc[nodes,"naphtha"].sum()/8760.)


    network.madd("Load",
                 ["kerosene for aviation"],
                 bus="EU oil",
                 carrier="kerosene for aviation",
                 p_set = get_parameter(options["aviation_demand"])*nodal_energy_totals.loc[nodes,["total international aviation","total domestic aviation"]].sum(axis=1).sum()*1e6/8760.)

    #NB: CO2 gets released again to atmosphere when plastics decay or kerosene is burned
    #except for the process emissions when naphtha is used for petrochemicals, which can be captured with other industry process emissions
    #tco2 per hour
    co2 = network.loads.loc[["naphtha for industry","kerosene for aviation","oil for shipping"],"p_set"].sum()*costs.at["oil",'CO2 intensity'] - industrial_demand.loc[nodes,"process emission from feedstock"].sum()/8760.

    network.madd("Load",
                 ["oil emissions"],
                 bus="co2 atmosphere",
                 carrier="oil emissions",
                 p_set=-co2)

    network.madd("Load",
                 nodes,
                 suffix=" low-temperature heat for industry",
                 bus=[node + " urban central heat" if node + " urban central heat" in network.buses.index else node + " services urban decentral heat" for node in nodes],
                 carrier="low-temperature heat for industry",
                 p_set=industrial_demand.loc[nodes,"low-temperature heat"]/8760.)

    #remove today's industrial electricity demand by scaling down total electricity demand
    for ct in n.buses.country.dropna().unique():
        loads = n.loads.index[(n.loads.index.str[:2] == ct) & (n.loads.carrier == "electricity")]
        if n.loads_t.p_set[loads].empty: continue
        factor = 1 - industrial_demand.loc[loads,"current electricity"].sum()/n.loads_t.p_set[loads].sum().sum()
        n.loads_t.p_set[loads] *= factor

    network.madd("Load",
                 nodes,
                 suffix=" industry electricity",
                 bus=nodes,
                 carrier="industry electricity",
                 p_set=industrial_demand.loc[nodes,"electricity"]/8760.)

    network.madd("Bus",
                 ["process emissions"],
                 location="EU",
                 carrier="process emissions")

    #this should be process emissions fossil+feedstock
    #then need load on atmosphere for feedstock emissions that are currently going to atmosphere via Link Fischer-Tropsch demand
    network.madd("Load",
                 ["process emissions"],
                 bus="process emissions",
                 carrier="process emissions",
                 p_set = -industrial_demand.loc[nodes,["process emission","process emission from feedstock"]].sum(axis=1).sum()/8760.)

    network.madd("Link",
                 ["process emissions"],
                 bus0="process emissions",
                 bus1="co2 atmosphere",
                 carrier="process emissions",
                 p_nom_extendable=True,
                 efficiency=1.)

    #assume enough local waste heat for CC
    network.madd("Link",
                 ["process emissions CC"],
                 bus0="process emissions",
                 bus1="co2 atmosphere",
                 bus2="co2 stored",
                 carrier="process emissions CC",
                 p_nom_extendable=True,
                 capital_cost=costs.at["cement capture","fixed"],
                 efficiency=(1-costs.at["cement capture","capture_rate"]),
                 efficiency2=costs.at["cement capture","capture_rate"],
                 lifetime=costs.at['cement capture','lifetime'])



def add_waste_heat(network):

    print("adding possibility to use industrial waste heat in district heating")

    #AC buses with district heating
    urban_central = network.buses.index[network.buses.carrier == "urban central heat"]
    if not urban_central.empty:
        urban_central = urban_central.str[:-len(" urban central heat")]

        if options['use_fischer_tropsch_waste_heat']:
            network.links.loc[urban_central + " Fischer-Tropsch","bus4"] = urban_central + " urban central heat"
            network.links.loc[urban_central + " Fischer-Tropsch","efficiency4"] = 0.95 - network.links.loc[urban_central + " Fischer-Tropsch","efficiency"]

        for o in opts:
            if "B" in o:
                if options['use_biofuel_waste_heat']:
                    network.links.loc[urban_central + " biomass to liquid","bus4"] = urban_central + " urban central heat"
                    network.links.loc[urban_central + " biomass to liquid","efficiency4"] = 0.95 - network.links.loc[urban_central + " biomass to liquid","efficiency"]
                    network.links.loc[urban_central + " solid biomass to gas","bus4"] = urban_central + " urban central heat"
                    network.links.loc[urban_central + " solid biomass to gas","efficiency4"] = 0.95 - network.links.loc[urban_central + " solid biomass to gas","efficiency"]

        if options['use_fuel_cell_waste_heat']:
            network.links.loc[urban_central + " H2 Fuel Cell","bus2"] = urban_central + " urban central heat"
            network.links.loc[urban_central + " H2 Fuel Cell","efficiency2"] = 0.95 - network.links.loc[urban_central + " H2 Fuel Cell","efficiency"]

def decentral(n):
    n.lines.drop(n.lines.index,inplace=True)
    n.links.drop(n.links.index[n.links.carrier.isin(["DC","B2B"])],inplace=True)

def remove_h2_network(n):

    nodes = pop_layout.index

    n.links.drop(n.links.index[n.links.carrier.isin(["H2 pipeline"])],inplace=True)

    n.stores.drop(["EU H2 Store"],inplace=True)

    if options['hydrogen_underground_storage']:
        h2_capital_cost = costs.at["gas storage","fixed"]
        #h2_capital_cost = costs.at["hydrogen underground storage","fixed"]
    else:
        h2_capital_cost = costs.at["hydrogen storage","fixed"]

    #put back nodal H2 storage
    n.madd("Store",
           nodes + " H2 Store",
           bus=nodes + " H2",
           e_nom_extendable=True,
           e_cyclic=True,
           carrier="H2 Store",
           capital_cost=h2_capital_cost)


def add_biomass_transport(network):

    # costs for biomass transport
    transport_costs = pd.read_csv(snakemake.input.biomass_transport,
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

    network.madd("Link",
           biomass_transport.index,
           bus0=biomass_transport.bus0 + " solid biomass",
           bus1=biomass_transport.bus1 + " solid biomass",
           p_nom_extendable=True,
           length=biomass_transport.length.values,
           marginal_cost=biomass_transport.costs*0.01 * biomass_transport.length.values,
           capital_cost=1,
           carrier="solid biomass transport")


def get_parameter(item):
    """Check whether it depends on investment year"""
    if type(item) is dict:
        return item[investment_year]
    else:
        return item


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
           efficiency=1-0.03*n.lines.length/1000,
           marginal_cost=0,
           carrier="DC",
           length=n.lines.length,
           capital_cost=n.lines.capital_cost)

    # Remove AC lines
    print("Removing AC lines")
    lines_rm = n.lines.index
    n.mremove("Line", lines_rm)

    # Set efficiency of all DC links to include losses depending on length
    n.links.loc[n.links.carrier == 'DC', 'efficiency'] = 1 - 0.03 * n.links.loc[n.links.carrier == 'DC', 'length'] / 1000


#%%
if __name__ == "__main__":
    # Detect running outside of snakemake and mock snakemake for testing
    if 'snakemake' not in globals():
        from vresutils.snakemake import MockSnakemake
        snakemake = MockSnakemake(
            wildcards=dict(network='elec', simpl='', clusters='37', lv='1.0',
                           opts='', planning_horizons='2020',
                           sector_opts='120H-T-H-B-I-onwind+p3-dist1-cb48be3'),

            input=dict( network='../pypsa-eur/networks/elec_s{simpl}_{clusters}_ec_lv{lv}_{opts}.nc',
                        energy_totals_name='resources/energy_totals.csv',
                        co2_totals_name='resources/co2_totals.csv',
                        transport_name='resources/transport_data.csv',
                	    traffic_data = "data/emobility/",
                        biomass_potentials='resources/biomass_potentials.csv',
                        timezone_mappings='data/timezone_mappings.csv',
                        heat_profile="data/heat_load_profile_BDEW.csv",
                        costs="../technology-data/outputs/costs_{planning_horizons}.csv",
                	    h2_cavern = "data/hydrogen_salt_cavern_potentials.csv",
                        profile_offwind_ac="../pypsa-eur/resources/profile_offwind-ac.nc",
                        profile_offwind_dc="../pypsa-eur/resources/profile_offwind-dc.nc",
                        busmap_s="../pypsa-eur/resources/busmap_elec_s{simpl}.csv",
                        busmap="../pypsa-eur/resources/busmap_elec_s{simpl}_{clusters}.csv",
                        clustered_pop_layout="resources/pop_layout_elec_s{simpl}_{clusters}.csv",
                        simplified_pop_layout="resources/pop_layout_elec_s{simpl}.csv",
                        industrial_demand="resources/industrial_energy_demand_elec_s{simpl}_{clusters}.csv",
                        heat_demand_urban="resources/heat_demand_urban_elec_s{simpl}_{clusters}.nc",
                        heat_demand_rural="resources/heat_demand_rural_elec_s{simpl}_{clusters}.nc",
                        heat_demand_total="resources/heat_demand_total_elec_s{simpl}_{clusters}.nc",
                        temp_soil_total="resources/temp_soil_total_elec_s{simpl}_{clusters}.nc",
                        temp_soil_rural="resources/temp_soil_rural_elec_s{simpl}_{clusters}.nc",
                        temp_soil_urban="resources/temp_soil_urban_elec_s{simpl}_{clusters}.nc",
                        temp_air_total="resources/temp_air_total_elec_s{simpl}_{clusters}.nc",
                        temp_air_rural="resources/temp_air_rural_elec_s{simpl}_{clusters}.nc",
                        temp_air_urban="resources/temp_air_urban_elec_s{simpl}_{clusters}.nc",
                        cop_soil_total="resources/cop_soil_total_elec_s{simpl}_{clusters}.nc",
                        cop_soil_rural="resources/cop_soil_rural_elec_s{simpl}_{clusters}.nc",
                        cop_soil_urban="resources/cop_soil_urban_elec_s{simpl}_{clusters}.nc",
                        cop_air_total="resources/cop_air_total_elec_s{simpl}_{clusters}.nc",
                        cop_air_rural="resources/cop_air_rural_elec_s{simpl}_{clusters}.nc",
                        cop_air_urban="resources/cop_air_urban_elec_s{simpl}_{clusters}.nc",
                        solar_thermal_total="resources/solar_thermal_total_elec_s{simpl}_{clusters}.nc",
                        solar_thermal_urban="resources/solar_thermal_urban_elec_s{simpl}_{clusters}.nc",
                        solar_thermal_rural="resources/solar_thermal_rural_elec_s{simpl}_{clusters}.nc",
                	    retro_cost_energy = "resources/retro_cost_elec_s{simpl}_{clusters}.csv",
                        floor_area = "resources/floor_area_elec_s{simpl}_{clusters}.csv"
            ),
            output=['results/version-cb48be3/prenetworks/elec_s{simpl}_{clusters}_lv{lv}__{sector_opts}_{planning_horizons}.nc']
        )
        import yaml
        with open('config.yaml', encoding='utf8') as f:
            snakemake.config = yaml.safe_load(f)


    logging.basicConfig(level=snakemake.config['logging_level'])

    timezone_mappings = pd.read_csv(snakemake.input.timezone_mappings,index_col=0,squeeze=True,header=None)

    options = snakemake.config["sector"]

    opts = snakemake.wildcards.sector_opts.split('-')

    print('Options: ', opts)

    investment_year=int(snakemake.wildcards.planning_horizons[-4:])

    n = pypsa.Network(snakemake.input.network,
                      override_component_attrs=override_component_attrs)

    Nyears = n.snapshot_weightings.sum()/8760.

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout,index_col=0)
    pop_layout["ct"] = pop_layout.index.str[:2]
    ct_total = pop_layout.total.groupby(pop_layout["ct"]).sum()
    pop_layout["ct_total"] = pop_layout["ct"].map(ct_total.get)
    pop_layout["fraction"] = pop_layout["total"]/pop_layout["ct_total"]

    simplified_pop_layout = pd.read_csv(snakemake.input.simplified_pop_layout,index_col=0)

    costs = prepare_costs(snakemake.input.costs,
                          snakemake.config['costs']['USD2013_to_EUR2013'],
                          snakemake.config['costs']['discountrate'],
                          Nyears,
                          snakemake.config['costs']['lifetime'])

    remove_elec_base_techs(n)

    n.loads["carrier"] = "electricity"

    remove_non_electric_buses(n)

    n.buses["location"] = n.buses.index

    update_wind_solar_costs(n, costs)

    if snakemake.config["foresight"]=='myopic':
        add_lifetime_wind_solar(n)
        add_carrier_buses(n,snakemake.config['existing_capacities']['conventional_carriers'])

    add_co2_tracking(n)

    add_generation(n)

    add_storage(n)

    for o in opts:
        if o[:4] == "wave":
            wave_cost_factor = float(o[4:].replace("p",".").replace("m","-"))
            print("Including wave generators with cost factor of", wave_cost_factor)
            add_wave(n, wave_cost_factor)
        if o[:4] == "dist":
            snakemake.config["sector"]['electricity_distribution_grid'] = True
            snakemake.config["sector"]['electricity_distribution_grid_cost_factor'] = float(o[4:].replace("p",".").replace("m","-"))
        # if o == "bioT":
        #     options["biomass_transport"] = True

    nodal_energy_totals, heat_demand, ashp_cop, gshp_cop, solar_thermal, transport, avail_profile, dsm_profile, co2_totals, nodal_transport_data = prepare_data(n)

    if "nodistrict" in opts:
        options["central"] = False

    if "H" in opts:
        add_heat(n)

    if "I" in opts:
        add_industry(n)

    if options["hvdc"]:
        hvdc_transport_model(n)

    for o in opts:
        if "B" in o:
            add_biomass(n)
            if options["bioT"]:
                add_biomass_transport(n)

    if "T" in opts:
        add_land_transport(n)

    if "I" in opts and "H" in opts:
        add_waste_heat(n)

    if options['dac']:
        add_dac(n)

    if "decentral" in opts:
        decentral(n)

    if "noH2network" in opts:
        remove_h2_network(n)

    for o in opts:
        m = re.match(r'^\d+h$', o, re.IGNORECASE)
        if m is not None:
            n = average_every_nhours(n, m.group(0))
            break
    else:
        logger.info("No resampling")

    #process CO2 limit
    limit = get_parameter(snakemake.config["co2_budget"])
    print("CO2 limit set to",limit)

    for o in opts:

        if "cb" in o:
            path_cb = snakemake.config['results_dir'] + snakemake.config['run'] + '/csvs/'
            if not os.path.exists(path_cb):
                os.makedirs(path_cb)
            try:
                CO2_CAP=pd.read_csv(path_cb + 'carbon_budget_distribution.csv', index_col=0)
            except:
                build_carbon_budget(o)
                CO2_CAP=pd.read_csv(path_cb + 'carbon_budget_distribution.csv', index_col=0)

            limit=CO2_CAP.loc[investment_year]
            print("overriding CO2 limit with scenario limit",limit)

    for o in opts:
        if "Co2L" in o:
            limit = o[o.find("Co2L") + 4:]
            limit = float(limit.replace("p", ".").replace("m", "-"))
            print("overriding CO2 limit with scenario limit", limit)


    print("adding CO2 budget limit as per unit of 1990 levels of",limit)
    add_co2limit(n, Nyears, limit)

    for o in opts:

        if o[:10] == 'linemaxext':
            maxext = float(o[10:])*1e3
            print("limiting new HVAC and HVDC extensions to",maxext,"MW")
            n.lines['s_nom_max'] = n.lines['s_nom'] + maxext
            hvdc = n.links.index[n.links.carrier == 'DC']
            n.links.loc[hvdc,'p_nom_max'] = n.links.loc[hvdc,'p_nom'] + maxext


    if snakemake.config["sector"]['electricity_distribution_grid']:
        insert_electricity_distribution_grid(n)



    for o in opts:
        if "+" in o:
            oo = o.split("+")
            carrier_list=np.hstack((n.generators.carrier.unique(), n.links.carrier.unique(),
                                    n.stores.carrier.unique(), n.storage_units.carrier.unique()))
            suptechs = map(lambda c: c.split("-", 2)[0], carrier_list)
            if oo[0].startswith(tuple(suptechs)):
                carrier = oo[0]
                attr_lookup = {"p": "p_nom_max", "c": "capital_cost"}
                attr = attr_lookup[oo[1][0]]
                factor = float(oo[1][1:])
                #beware if factor is 0 and p_nom_max is np.inf, 0*np.inf is nan
                if carrier == "AC":  # lines do not have carrier
                    n.lines[attr] *= factor
                else:
                    comps = {"Generator", "Link", "StorageUnit"} if attr=='p_nom_max' else {"Generator", "Link", "StorageUnit", "Store"}
                    for c in n.iterate_components(comps):
                        if carrier=='solar':
                            sel = c.df.carrier.str.contains(carrier) & ~c.df.carrier.str.contains("solar rooftop")
                        else:
                            sel = c.df.carrier.str.contains(carrier)
                        c.df.loc[sel,attr] *= factor
                print("changing", attr ,"for",carrier,"by factor",factor)


    if snakemake.config["sector"]['gas_distribution_grid']:
        insert_gas_distribution_costs(n)
    if snakemake.config["sector"]['electricity_grid_connection']:
        add_electricity_grid_connection(n)

    n.export_to_netcdf(snakemake.output[0])
