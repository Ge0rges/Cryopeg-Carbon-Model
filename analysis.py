import numpy as np
from plots import plot_model, plot_sensitivity
from julia import Main
from decimal import Decimal
from sympy import symbols, exp, integrate


Main.include("model.jl")


## ASSUMPTIONS
# Assumed each cell has equal carbon content.
# Assumed organic carbon content is equal throughout cell lifespan, including death.
# Assumed cells contain no inorganic carbon upon death.
# Assumed constant death fraction.
# Assumed each live cell divides at same rate.
# Assumed constant carbon fixing factor.
# Assumed organic carbon content arround the brine is equal throughout brine lifespan.
# Assumed diffusion has no effect within the brine.
# Assumed no external disturbance to the brine.
# Assumed selected organism represents entire brine microbial community kinetics.
# Assumed organic carbon to be the only limiting element.
## -> Therefore, assumed the carbon source is also providing nitrogen.
# Assumed dry mass of cell is equal to carbon content per cell.
# Assumed that only organic carbon towards maintenance is respired.
# Assumed no other organic carbon costs than maintenance and dndt * resource_quantity_per_individual
######
# Assumed brine is a prolate spheroid.
# Assumes proportional expansion of the brine volume in every direction.
# Assumes brine volume expansion cycle occurs once a year, in one peak.


## MODEL PARAMATERS
default_paramaters = [10 ** 9, 0, 0, 0.001, 0, 882000000, 140, 0, 0.06, 0.0006/(0.00069 * 0.4), 0]
default_paramater_names = ["carrying_capacity", "organic_carbon_input", "inorganic_carbon_input",
                           "Death %", "IC F%", "Ks", "OC/Cell", "IC/Cell", "mu_max", "ME", "m_prime"]
default_paramater_bounds = [None, None, None, None, [0, 1], [0.001, 10**10], [0, 500], [0, 500], [0.000000001, 0.1], [0, 400], None]
assert len(default_paramaters) == len(default_paramater_names) == len(default_paramater_bounds)

# carrying_capacity inncells/ml - Maximum cell density.
# organic_carbon_input in fg/(ml x day) - Organic arbon input per day.
# inorganic_carbon_input in "fg/(ml x day) - Inorganic carbon input per day.
# natural_death_fraction in %/day - Fraction of cells present that will die at any time t.
# inorganic_carbon_fixing_factor in %/day - The percentage of inorganic carbon present fixed
# ks in fg C - The organic carbon concentration at which u = 1/2 u0.
# organic_carbon_content_per_cell in fg/cell -  Amount of organic carbon per cell.
# inorganic_carbon_content_per_cell in fg/cell - Amount of inorganic carbon per cell.
# mu_max in day^-1 - Max growth rate.
# base_maintenance_per_cell in fg carbon/fg dry mass x day - Constant maintenance energy coefficient
# m_prime in fg glucose/fg dry mass x day - Growth rate dependent maintenance coefficient, m_prime in pirt (1982)

default_ivp = [1286820000000, 0, 100000, 11000 * 365.25]
default_ivp_names = ["initial_carbon_content", "initial_inorganic_carbon_content", "initial_cell_count", "duration"]

# initial_carbon_content in fg/ml - Initial amount of organic carbon in the system per ml
# initial_inorganic_carbon_content in fg/ml - Initial amount of inorganic carbon in the system per ml.
# initial_cell_count in cells/ml - Initial number of cells in the system per ml.
# duration in days - Duration of the system.


def estimate_me_no_growth(start_carbon, end_carbon, end_cell, timespan):
    # start_cell in cells/ml - Start cell abundance per ml
    # end_cell in cells/ml - End cell abundance per ml
    # timespan in days - Community lifespan
    # start_carbon in fg C/ml - Start carbon concentration
    # end_carbon in fg C/ml - End carbon concentration

    m = (end_carbon - start_carbon) / (end_cell * timespan)
    return m


def estimate_me_exp_growth(start_carbon, end_carbon, start_cell, end_cell, cell_carbon_content, timespan):
    # start_cell in cells/ml - Start cell abundance per ml
    # end_cell in cells/ml - End cell abundance per ml
    # timespan in days - Community lifespan
    # start_carbon in fg C/ml - Start carbon concentration
    # end_carbon in fg C/ml - End carbon concentration
    # cell_carbon_content in fg C - Carbon content a cell takes up/releases upon birth/death

    # N(t) = N_0 * exp(μ*t) where μ is growth rate
    N0, mu, t = symbols("N_0 mu t", real=True)

    N0 = start_cell
    mu = np.log(end_cell / start_cell) / timespan  # days ^-1
    # number_of_generations = np.log(end_cell / start_cell) / np.log(2)  # generations
    # doubling_time = timespan/number_of_generations  # days

    N = N0 * exp(mu * t)

    # dC/dt = -g*dN/dt - Nm   where: g is cell_carbon_content,
    # m is maintenance energy, C is carbon concentration, N is cell concentration
    # Integrate from 0 to timespan.
    # Cf - C0 = -g(Nf-N0) - m∫N(t)dt
    # m = (g(Nf-N0) + Cf - C0)/∫N(t)dt
    m = (cell_carbon_content * (end_cell - start_cell) + end_carbon - start_carbon) / integrate(N, (t, 0, timespan))

    return mu, m


def run_model(p, ivp):
    # Run the model trhough PyJulia
    S, I, N, t = Main.run_model(default_paramaters, default_ivp)
    return S, I, N, t


def run_sensitivity_analysis(p, p_bounds, ivp):
    # Sensitivity analysis
    # Fix the bounds to replace "None" with no range bounds tha twill yield a 0 sobol index.
    for i, x in enumerate(p_bounds):
        if x is None:
            p_bounds[i] = [p_bounds[i], p_bounds[i]]

    # Run the sensitivity analysis
    ST, S1 = Main.run_sensitivity_analysis(p, p_bounds, ivp)
    return ST, S1


def calculate_brine_expansion(carbon_density_in_permafrost, carbon_required_per_year):
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
    ratio = volume_needed / brine_volume

    # Expansion ratio of area
    ratio_area = (1 + volume_needed / brine_volume) ** Decimal('2')/Decimal('3')

    # Expansion ratio linearly
    ratio_dimensions = (1 + volume_needed / brine_volume) ** Decimal('1')/Decimal('3')

    expansion_a = a * ratio_dimensions
    expansion_b = b * ratio_dimensions

    return ratio, ratio_area, ratio_dimensions, expansion_a, expansion_b


# Example code for all analysis on default values
if __name__ == "__main__":
    # Estimate bounds for maintenance energy
    me_upper = estimate_me_no_growth(default_ivp[1], default_ivp[1]/10, 10 ** 8, default_ivp[3])
    growth_rate, me_lower = estimate_me_exp_growth(default_ivp[1], default_ivp[1]/10, 10**5, 10**8, 400, default_ivp[3])
    print("Bounds for maintenance energy are, low: " + str(me_lower) + " and high: " + str(me_upper) + "fg C/cell day")

    # Run the model trhough PyJulia
    S, I, N, t = run_model(default_paramaters, default_ivp)
    plot_model(S, I, N, t)

    # Calculate brine expansion
    _, _, ratio_dimensions, _, _ = calculate_brine_expansion(default_ivp[0], (S[-1] - S[0])/(default_ivp[3] * 365.25))
    print("Brine needs to expand by " + str(ratio_dimensions) + "% along each axis per year.")

    # Run the sensitivity analysis
    ST, S1 = Main.run_sensitivity_analysis(default_paramaters, default_paramater_bounds, default_ivp)
    plot_sensitivity(ST, S1, default_paramater_names)
