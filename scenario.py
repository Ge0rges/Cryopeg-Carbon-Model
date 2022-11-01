"""
Defines the Scenario class.
Functions which instantiate CB1, CB4 and CBIW scenarios.
"""

from utils import *


class Scenario:
    """
    Scenarios encapsulate model constants, variable inputs, and model results.
    """

    scenario_name = None  # string - Used to title plots resulting from this scenario.

    # Variable paramaters
    start_poc = None  # fg C - The measured particulate organic carbon around the brine today.
    end_poc = None  # fg C - The measured particulate organic carbon in the brine today.
    start_doc = None  # fg C - The measured dissolved organic carbon around the brine today.
    end_doc = None  # fg C - The measured dissolved organic carbon in the brine today.

    punctual_organic_carbon_addition = []  # (time in days, (fg pOC/ml, fg dOC/ml)) - Added carbon at a given time.

    observed_end_cell_density = None  # cells - The measured cell density in the corresponding brine today.

    lab_growth_rate = None  # 1/day - Growth rate determined under nutrient replete in-situ conditions in the lab.
    maintenance_per_cell = None  # fg C/day cell - Constant maintenance energy coefficient.
    dissolved_organic_carbon_per_cell = None  # fg dOC/cell -  Amount of dissolved organic carbon per cell.

    title = None  # String - Scenario name. Used in plots later.

    # Constants - Common to all Scenarios
    _start_cell = 10 ** 5  # Average cells/ml order of magnitue for sea-ice (Cooper et al. 2019)
    _carrying_capacity = 10 ** 9  # cells/ml - Maximum cell density.
    _growth_rate = None  # 1/day - Growth rate used by the model.

    _particulate_organic_carbon_input_rate = 0  # fg C/(ml x day) - Particulate organic carbon input per day.
    _dissolved_organic_carbon_input_rate = 0  # fg C/(ml x day) - Dissolved carbon input per day.

    _start_inorganic_carbon_content = 0  # fg C/ml - Inroagnic carbon in brine at start.
    _inorganic_carbon_input_rate = 0  # fg C/(ml x day) - Organic carbon input per day.
    _inorganic_carbon_fixation_rate = 0  # %/day - Percentage of inorganic carbon converted to organic carbon per day.
    _inorganic_carbon_per_cell = 0  # fg C/cell -  Amount of inorganic carbon per cell.

    _timespan = 40000*365.25  # Age in days based on carbon dating (Iwanaha et al. 2021)

    _eea_rate = 0.012192  # fg pOC/cell * day - The hydrolysis rate of pOC by extracellular enzymes (Showalter, 2021)

    # _Kd: fg C - The organic carbon concentration at which u = 1/2 u0.
    _Kd = 8.82 * 10 ** 5  # average Ks ug AA/L = 2.152 (Yager & Deming 1999) * DCAA are 41% carbon (Rowe & Deming 1985)

    _paramater_bounds = [[1e-6, 1e2], [1e-5, 500], [1, 500], None, None,  # Ordered bounds for sensitivty analysis
                         None, None, None, None, [0, 100], [1e3, 1e8], None]  # as in julia p list.

    _paramater_names = ["μ", "Maintenance energy", "dOC/cell", "Carrying capacity", "pOC input rate",
                        "dOC input rate", "Inorganic carbon input rate", "Inorganic carbon fixation rate", "IC/cell",
                        "EEA rate", "Kd", "Punctual organic carbon addition"]

    # Methods
    def get_julia_ordered_paramaters(self):
        """
        The julia model takes in the paramaters as an ordered list. This method builds and returns that list.
        :return: Ordered parameter list
        :rtype: List
        """
        punctual_organic_carbon_addition = -1 if len(self.punctual_organic_carbon_addition) == 0 else self.punctual_organic_carbon_addition  # Can't pass empty arrays to Julia

        ordered_p = [self._growth_rate, self.maintenance_per_cell, self.dissolved_organic_carbon_per_cell,
                     self._carrying_capacity, self._particulate_organic_carbon_input_rate,
                     self._dissolved_organic_carbon_input_rate, self._inorganic_carbon_input_rate,
                     self._inorganic_carbon_fixation_rate, self._inorganic_carbon_per_cell,
                     self._eea_rate, self._Kd, punctual_organic_carbon_addition]
        # punctual_organic_carbon_addition always at end in julia code. Changes to this array must be manually carried
        # to model.jl run_model()

        assert len(self._paramater_bounds) == len(self._paramater_names) == len(ordered_p)

        return ordered_p

    def get_julia_ordered_ivp(self):
        """
        The julia model takes in the initial variable problem as an ordered list.
        This method builds and returns that list.
        :return: Ordered ivp list
        :rtype: List
        """

        return [self.start_poc, self.start_doc, self._start_inorganic_carbon_content, self._start_cell, self._timespan]


def cb1_scenario():
    """
    CB1 is an intra-sediment brine for which specific carbon values of its immediate surroundings are unavailable.
    However, it shares many characteristics with the other intra-sediment brines studied. Therefore, using its data,
    paired with regional organic carbon values, we construct a common intra-sediment brine scenario.
    """
    dry_sediment_density = 2.625 * 10 ** 6  # ug/ml average density density of kaolinite and sand
    expansion_factor = 0.0905  # Assumed expansion factor of porewater from liquid to solid. (French & Shur 2010)
    beo_volumetric_ice_content = 0.731  # L solid porewater/L permafrost from Go Iwahana (unpublished, but see Meyer et al. 2010)
    beo_poc = 0.0232  # ug C/ug dry sed (our data)
    beo_doc = 51.2323  # mg C/L thawed porewater (our data)
    cb1_18_poc = 1.24 * 10 ** 4  # in uM from Cooper et al. 2019
    cb1_18_doc = 1.02 * 10 ** 5  # in uM from Cooper et al. 2019

    # Conversions
    # Converting from ug C/ug dry sediment to fg C/ml permafrost
    # L dry sediment / L permafrost = (1 - L solid porewater / L permafrost)
    # ug C/ug dry sed * ug dry sediment/ml dry sediment * ml/dry sediment/ml permafrost
    beo_poc *= dry_sediment_density * (1 - beo_volumetric_ice_content) * 10 ** 9

    # Converting from mg C/L thawed porewater * L thawed porewater/L permafrost to fg C/ml
    beo_doc *= beo_volumetric_ice_content * (1 - expansion_factor) * 10 ** 9

    # Converting from micromolar C to fg C/ml using carbon molar mass of 12.011
    cb1_18_poc = convert_micromolar_carbon_to_fgC_per_ml(cb1_18_poc)
    cb1_18_doc = convert_micromolar_carbon_to_fgC_per_ml(cb1_18_doc)

    # PARAMATERS
    scenario = Scenario()
    scenario.title = "CB1"

    scenario.start_poc = beo_poc
    scenario.start_doc = beo_doc
    scenario.end_poc = cb1_18_poc
    scenario.end_doc = cb1_18_doc
    scenario.observed_end_cell_density = 5.7 * 10 ** 6  # Average cell/ml order of magnitude (Cooper et al. 2019)
    scenario.lab_growth_rate = 0.06  # Marinobacter aerobic growth rate in-situ based on lab experiments (unpublished), 50% anaerobic penalty
    scenario.dissolved_organic_carbon_per_cell = 15.7  # Litterature based value (Nguyen & Maranger 2011)

    return scenario


def cb4_scenario():
    """
    CB4 is an intra-sediment brine for which carbon values of brine and its immediate surrounding are available. This
    is therefore a specific case of the more general "intra-sediment" brine scenario. Its biology is also somewhat
    unique from the other cryopeg brines, e.g. Psychrobacter is the dominant genus in this brine rather than
    Marinobacter like many others.
    """
    dry_sediment_density = 2.625 * 10 ** 6  # ug/ml average density density of kaolinite and sand
    expansion_factor = 0.0905  # Assumed expansion factor of porewater from liquid to solid. (French & Shur 2010)
    cb4_volumetric_ice_content = 0.527  # L solid porewater/L permafrost from Go Iwahana (unpublished, but see Meyer et al. 2010)
    cb4_poc = 0.0136  # ug C/ug dry sed (our data)
    cb4_doc = 1286.81  # mg C/L thawed porewater (our data)
    cb4_brine_poc = 4.14 * 10 ** 3  # in uM (Cooper et al. 2019)
    cb4_brine_doc = 8.50 * 10 ** 4  # in uM (Cooper et al. 2019)

    # Conversions
    # Converting from ug C/ug dry sediment to fg C/ml permafrost
    # L dry sediment / L permafrost = (1 - L solid porewater / L permafrost)
    # ug C/ug dry sed * ug dry sediment/ml dry sediment * ml/dry sediment/ml permafrost
    cb4_poc *= dry_sediment_density * (1 - cb4_volumetric_ice_content) * 10 ** 9

    # Converting from mg C/L thawed porewater * L thawed porewater/L permafrost  to fg C/ml
    cb4_doc *= cb4_volumetric_ice_content * (1 - expansion_factor) * 10 ** 9

    # Converting from micromolar C to fg C/ml using carbon molar mass of 12.011
    cb4_brine_poc = convert_micromolar_carbon_to_fgC_per_ml(cb4_brine_poc)
    cb4_brine_doc = convert_micromolar_carbon_to_fgC_per_ml(cb4_brine_doc)

    # PARAMATERS
    scenario = Scenario()
    scenario.title = "CB4"

    scenario.start_poc = cb4_poc
    scenario.start_doc = cb4_doc
    scenario.end_poc = cb4_brine_poc
    scenario.end_doc = cb4_brine_doc
    scenario.observed_end_cell_density = 1.14 * 10 ** 7  # cell/ml for CB4_18 (Cooper et al. 2019)
    scenario.lab_growth_rate = 0.016  # Psychrobacter cryohalolentis growth rate in-situ based on lab experiments (Bakermans et al. 2003)
    scenario.dissolved_organic_carbon_per_cell = 54.04  # Took average P. cryohalolentis size: 0.365014 um3 (Bakermans et al. 2006), and carbon conversion factor as 148 fg C/um3 (Kirchman et al. 2009)

    return scenario


def cbiw_scenario():
    """
    CBIW is a intra-ice cryopeg brine. It is thought to have migrated from sediment to massive ice at around 11000 years
    ago. Therefore, this scenario is unique in including an injection of carbon at 11 000 years. We use the lower bound
    estimate of how much carbon was injected to be the organic carbon measured in the massive ice around the brine today.
    It is likely that more carbon was added than our estimate, however given our current understanding it is impossible
    to quantify precisely how much carbon was added into the brine during its migration.
    """
    dry_sediment_density = 2.625 * 10 ** 6  # ug/ml average density density of kaolinite and sand
    expansion_factor_porewater = 0.0905  # Assumed expansion factor of porewater from liquid to solid. (French & Shur 2010)
    beo_volumetric_ice_content = 0.731  # L solid porewater/L permafrost from Go Iwahana (unpublished, but see Meyer et al. 2010)
    expansion_factor_ice = 0.08042  # % Expansion factor of pure water (massive ice) from liquid to solid
    beo_poc = 0.0232  # ug C/ug dry sed (our data)
    beo_doc = 51.2323  # mg C/L thawed porewater (our data)
    massive_ice_poc = 4.5  # ug C/ml thawed massive ice (Collangelo-Lillis 2016)
    massive_ice_doc = 0.65  # μg C/ml thawed massive ice. EPS as a proxy for DOC. (Colangelo-Lillis et al. 2016)
    cbiw_brine_poc = 1.98 * 10 ** 3  # in uM (Cooper et al. 2019)
    cbiw_brine_doc = 3 * 10 ** 4  # in uM (Cooper et al. 2019)

    # Conversions
    # Converting from ug C/ug dry sediment to fg C/ml permafrost
    # L dry sediment / L permafrost = (1 - L solid porewater / L permafrost)
    # ug C/ug dry sed * ug dry sediment/ml dry sediment * ml/dry sediment/ml permafrost
    beo_poc *= dry_sediment_density * (1 - beo_volumetric_ice_content) * 10 ** 9

    # Converting from mg C/L thawed porewater * L thawed porewater/L permafrost to fg C/ml
    beo_doc *= beo_volumetric_ice_content * (1 - expansion_factor_porewater) * 10 ** 9

    # Converting from ug C/ml thawed massive ice to fg C/ml massive ice
    massive_ice_poc *= massive_ice_poc * (1 - expansion_factor_ice) * 10 ** 9
    massive_ice_doc *= massive_ice_doc * (1 - expansion_factor_ice) * 10 ** 9

    # Converting from micromolar C to fg C/ml using carbon molar mass of 12.011
    cbiw_brine_poc = convert_micromolar_carbon_to_fgC_per_ml(cbiw_brine_poc)
    cbiw_brine_doc = convert_micromolar_carbon_to_fgC_per_ml(cbiw_brine_doc)

    # Build the scenario
    scenario = Scenario()
    scenario.title = "CBIW"

    scenario.start_poc = beo_poc
    scenario.start_doc = beo_doc
    scenario.end_poc = cbiw_brine_poc
    scenario.end_doc = cbiw_brine_doc

    scenario.punctual_organic_carbon_addition = [
        (scenario._timespan - 11000*365.25, (massive_ice_poc, massive_ice_doc))  # Add massive ice carbon at 11000y ago.
    ]

    scenario.observed_end_cell_density = 1.39 * 10 ** 8  # average cell/ml of CBIW (Cooper et al. 2019)
    scenario.lab_growth_rate = 0.06  # Marinobacter aerobic growth rate in-situ based on lab experiments (unpublished), 50% anaerobic penalty
    scenario.dissolved_organic_carbon_per_cell = 15.7  # Litterature based value (Nguyen & Maranger 2011)

    return scenario
