"""
Contains functions that encapsulate or execute various analyses or calculations.
Functions take a scenario object, an
"""

import numpy as np
import copy

from utils import *
from scenario import Scenario
from julia import Main
from decimal import Decimal
from sympy import symbols, exp, integrate, log, Eq, solveset, Reals


Main.include("model.jl")


### MODEL ASSUMPTIONS ###
# Assumed each cell has equal carbon content.
# Assumed cell organic carbon content is equal throughout cell lifespan, including death.
# Assumed cells contain no inorganic carbon upon death.
# Assumed equal growth rate across the community.
# Assumed constant carbon fixing factor.
# Assumed organic carbon content arround the brine is equal throughout brine lifespan.
# Assumed instant diffusion of material across the brine.
# Assumed selected model organism represents entire brine microbial community kinetics.
# Assumed organic carbon to be the only limiting element.

### BRINE EXPANSION ASSUMPTIONS ###
# Assumed brine is a prolate spheroid.
# Assumes proportional expansion of the brine volume in every direction.
# Assumes brine volume expansion cycle occurs once a year, in one peak.


class Analysis:
    """
    A class that encapsulates analysis settings and results.
    """
    scenario: Scenario = None

    _use_minimum_growth_rate: bool = None
    _use_me_lower_bound: bool = None

    maintenance_energy_result: MaintenanceResult = None
    sensitivity_analysis_result: SensitivityResult = None
    model_result: ModelResult = None
    expansion_result: ExpansionResult = None
    end_eps: float = None
    eea_estimation: EEAResult = None

    title = None

    def __init__(self, scenario, use_minimum_growth_rate, use_me_lower_bound):
        """
        Initiates the analysis and construts its title.
        """
        self.scenario = copy.deepcopy(scenario)
        self._use_minimum_growth_rate = use_minimum_growth_rate
        self._use_me_lower_bound = use_me_lower_bound

        # Invalidate any growth rate present, we need to calcualte it
        if self._use_minimum_growth_rate:
            self.scenario.growth_rate = None

        self.title = self.scenario.title + " - "
        self.title += "minimum growth rate - " if use_minimum_growth_rate else "lab growth rate - "
        self.title += "lower bound ME" if use_me_lower_bound else "upper bound ME"

        return


    def run_analysis(self, do_sensitivity_analysis):
        """
        Runs all analyses including sensitivity analysis if requested. Configures the scenario's growth rate and
        maintenance energy as required by the analysis. Stores all results in the class.
        """
        # Estimate bounds for maintenance energy
        self.maintenance_energy_result = estimate_me_bounds(self.scenario)

        # Switch out growth rate to calcualted one if required
        self.scenario._growth_rate = self.maintenance_energy_result.minimum_growth_rate if self._use_minimum_growth_rate else self.scenario.lab_growth_rate

        # Switch out maintenance energy based on paramater
        self.scenario.maintenance_per_cell = self.maintenance_energy_result.lower_bound_me if self._use_me_lower_bound else self.maintenance_energy_result.upper_bound_me

        # Run the model through PyJulia using the lower ME
        self.model_result = run_model(self.scenario)

        # Calculate brine expansion
        self.expansion_result = calculate_brine_expansion(self.scenario.start_poc+self.scenario.start_doc,
                                                          (self.model_result.dOC[-1] - self.model_result.dOC[0]) / self.scenario._timespan)

        # Run the sensitivity analysis
        self.sensitivity_analysis_result = run_sensitivity_analysis(self.scenario) if do_sensitivity_analysis else None

        # Calculate eps produced
        self.end_eps = estimate_end_eps(self)

        # Calculate EEA rate bounds
        self.eea_estimation = estimate_eea_rate(self.scenario)


def estimate_me_bounds(scenario: Scenario):
    """
    Calculates the maintenance energy bounds and minimum growth rate of the system.
    This is done by assuming that organic carbon consumed is equal to the difference between total organic carbon
    in the brine, and (total organic brine surrounding the brine + any carbon additions).
    Upper bound corresponds to a strict expoential growth case from assumed start cell density to observed end cell
    density.
    Lower bound corrsponds to a case where no growth occured, meaning start cell density is equal to end cell density.
    Returns a tuple: (minimum growth rate, minimum doubling time, ME upper bound, ME lower bound).
    ME values are in fg C/cell day.
    """

    added_carbon_total = scenario._timespan * (scenario._particulate_organic_carbon_input_rate
                                              + scenario._dissolved_organic_carbon_input_rate)
    for p_add in scenario.punctual_organic_carbon_addition:
        added_carbon_total += p_add[1][0] + p_add[1][1]

    start_cell = scenario._start_cell
    end_cell = scenario.observed_end_cell_density
    timespan = scenario._timespan
    start_carbon = scenario.start_poc + scenario.start_doc + added_carbon_total
    end_carbon = scenario.end_poc + scenario.end_doc
    cell_carbon_content = scenario.dissolved_organic_carbon_per_cell

    # log(N_f/N_0) = μ*(tf-t0) where μ is growth rate
    N0, mu, t = symbols("N_0 mu t", real=True)  # t is time elapsed

    number_of_generations = np.log(end_cell / start_cell) / np.log(2)  # generations
    doubling_time = timespan / number_of_generations  # days
    mu = np.log(end_cell / start_cell) / timespan  # /days

    N0 = start_cell
    N = exp(mu * t + log(N0))

    # dS/dt = -g*dN/dt - Nm   where: g is cell_carbon_content,
    # m is maintenance energy, S is carbon concentration, N is cell concentration
    # Integrate from 0 to timespan.
    # Sf - S0 = -g(Nf-N0) - m∫N(t)dt
    # m = (S0 - Sf - g(Nf-N0))/∫N(t)dt
    m_up = (start_carbon - end_carbon - cell_carbon_content * (end_cell - start_cell)) / integrate(N, (t, 0, timespan))
    m_low = (start_carbon - end_carbon) / (end_cell * timespan)

    # Pack result
    result = MaintenanceResult()
    result.minimum_growth_rate = mu
    result.minimum_doubling_time = doubling_time
    result.lower_bound_me = m_low
    result.upper_bound_me = m_up

    return result


def estimate_end_eps(analysis: Analysis):
    """
    Estimates total EPS produced by the system given a model's cell over time result.
    """
    assert analysis.model_result

    eps_rate = 0.337712  # fg C/cell day From Stoderegger and Herndl, 1998 (7.1 amol C/cell hr). Applied Q10 (f=2), 20C.

    cells = analysis.model_result.cells
    timepoints = analysis.model_result.t

    total_eps = 0

    for i, cell_count in enumerate(cells):
        if i < len(cells) - 1 and cell_count >= 1:
            delta_t = timepoints[i+1] - timepoints[i]
            total_eps += delta_t*eps_rate*cell_count

    return total_eps


def estimate_eea_rate(scenario: Scenario):
    """
    Estimates the extracellular enzyme activity rate for the total amount of POC converted to DOC in a given timeframe.
    """

    start_cell = scenario._start_cell
    end_cell = scenario.observed_end_cell_density
    timespan = scenario._timespan

    # Get total amount of POC added
    total_poc_added = scenario._timespan * scenario._particulate_organic_carbon_input_rate

    for p_add in scenario.punctual_organic_carbon_addition:
        total_poc_added += p_add[1][0]

    # Get total amount of POC converted to DOC
    poc_converted = (scenario.start_poc + total_poc_added) - scenario.end_poc

    # Build growth curve for integral
    # log(N_f/N_0) = μ*(tf-t0) where μ is growth rate
    N0, mu, t, p = symbols("N_0 mu t p", real=True)  # t is time elapsed

    mu = np.log(scenario.observed_end_cell_density / start_cell) / timespan  # /days

    N0 = start_cell
    N = exp(mu * t + log(N0))  # Cell growth function

    # Calculate EEA bounds considering two growth cases, as in maintenance energy.
    eea_upper = poc_converted/integrate(N, (t, 0, timespan))
    eea_lower = poc_converted/(end_cell * timespan)

    # Calculate the timespan the EEA rate in-use
    eq = Eq(integrate(N, (t, 0, p)), poc_converted/scenario._eea_rate)

    predicted_timespan = solveset(eq, p, domain=Reals)

    result = EEAResult()
    result.eea_upper = eea_upper
    result.eea_lower = eea_lower
    result.predicted_timespan = predicted_timespan.args[0]

    return result


def run_model(scenario: Scenario):
    """
    Runs the model by interfacing with the Julia code through PyCall. Model is run on passed scenario.
    Returns model output as tuple of lists: pOC, dOC, inorganic carbon, cell count, time
    """
    # Run the model through PyJulia
    P, D, I, N, t = Main.run_model(scenario.get_julia_ordered_paramaters(), scenario.get_julia_ordered_ivp())

    result = ModelResult()
    result.pOC = P
    result.dOC = D
    result.IC = I
    result.cells = N
    result.t = t

    return result


def run_sensitivity_analysis(scenario: Scenario):
    """
    Runs a sensitivity analysis of the model by interfacing with the Julia code through PyCall.
    Sensitivty analysis is run using passed Scenario. Note that variables with non-None bounds will be varied.
    Returns total and first order sobol indices in a tuple ordered as such.
    """

    p_bounds = scenario._paramater_bounds
    p = scenario.get_julia_ordered_paramaters()

    # Sensitivity analysis
    # Fix the bounds to replace "None" with no range bounds tha twill yield a 0 sobol index.
    for i, x in enumerate(p_bounds):
        if x is None:
            p_bounds[i] = [p[i], p[i]]

    # Run the sensitivity analysis
    ST, S1 = Main.run_sensitivity_analysis(p_bounds, scenario.get_julia_ordered_ivp(), False)

    result = SensitivityResult()
    result.total_sobol_indices = ST
    result.first_order_sobol_indices = S1

    return result


def calculate_brine_expansion(carbon_density_in_permafrost, carbon_required_per_year):
    """
    Given an organic carbon need per year, and the organic carbon density of permafrost surrounding a brine, this
    function calculates by what distance the brine (a prolate spheroid) would have to expand to incorporate that carbon.
    This assumes that each year the permafrost the brine was "thawed into" previously is replenished in carbon.
    This assumes a proportional expansion in all directions of the brine.
    """
    ## BRINE VOLUME PARAMATERS
    # Prolate spheroid dimensions
    a = Decimal('2.185099')  # cm - width
    b = Decimal('4')  # cm - length
    brine_volume = Decimal('4') / Decimal('3') * Decimal('3.14159') * (a / Decimal('2')) ** Decimal('2') * (
            b / Decimal('2'))  # Volume of prolate spheroid: V = 4/3 * pi * (a/2) ** 2 * (b/2)

    # Convert to Decimal objects for accuracy
    carbon_required_per_year = Decimal(carbon_required_per_year)
    carbon_density_in_permafrost = Decimal(carbon_density_in_permafrost)

    # Calculate volume needed to expand
    # The volume that must be covered to get the required amount of C assuming all carbon is replenished each year,
    # and carbon from permafrost "expanded into" is added to the brine with no loss.
    volume_needed = carbon_required_per_year / carbon_density_in_permafrost

    # Expansion ratio of volume
    ratio_volume = volume_needed / brine_volume

    # Expansion ratio of area
    ratio_area = (1 + volume_needed / brine_volume) ** Decimal('2') / Decimal('3')

    # Expansion ratio linearly
    ratio_dimensions = (1 + volume_needed / brine_volume) ** Decimal('1') / Decimal('3')

    expansion_a = a * ratio_dimensions
    expansion_b = b * ratio_dimensions

    # Pack result
    result = ExpansionResult()
    result.ratio_volume = ratio_volume
    result.ratio_area = ratio_area
    result.ratio_dimensions = ratio_dimensions
    result.expansion_a = expansion_a
    result.expansion_b = expansion_b

    return result
