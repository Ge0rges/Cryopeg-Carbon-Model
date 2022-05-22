from decimal import Decimal
from typing import List

import matplotlib.pyplot as plt
import numpy as np
from decimal import *
from scipy.integrate import solve_ivp
from sympy import symbols, Max, Eq, Function, classify_ode, diff
from sympy.solvers.ode.systems import dsolve_system

from SALib.sample import saltelli
from SALib.analyze import sobol


plt.rcParams['figure.figsize'] = 12, 10

getcontext().rounding = ROUND_CEILING
np.seterr(all='raise')

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


## CONSTANTS
# conversions
one_kilo_in_femto = Decimal('1e18')
one_microgram_in_femto = Decimal('1e9')
one_m2_in_cm2 = Decimal('10000')
one_gram_in_femto = Decimal('1e15')
one_gram_in_microgram = Decimal('1e6')

## BRINE VOLUME PARAMATERS
# Prolate spheroid dimensions
a = Decimal('2.185099')  # cm - width
b = Decimal('4')  # cm - length
brine_volume = Decimal('4') / Decimal('3') * Decimal('3.14159') * (a / Decimal('2')) ** Decimal('2') * (b / Decimal('2'))  # Volume of prolate spheroid: V = 4/3 * pi * (a/2) ** 2 * (b/2)

## MODEL PARAMATERS
# # Organic carbon around brine
# utgiagvik_nutrients_per_cm2_for_100cm = Decimal('83.4') / one_m2_in_cm2 * one_kilo_in_femto  # fg/cm2 - Initial C available per cm2 over 100cm depth if SOCC100.
# utgiagvik_nutrients_depth = 100  # cm - The NCSCD gives cm2 numbers over a total depth.

default_paramaters = {
    # Inputs
    "organic_carbon_input": [Decimal('0'), "fg/(ml x day) - Organic arbon input per day.", [0, 10**10]],
    "inorganic_carbon_input": [Decimal('0'), "fg/(ml x day) - Inorganic carbon input per day.", [0, 10**10]],

    # Community paramaters
    "natural_death_fraction": [Decimal('0.001'), "%/day - Fraction of cells present that will die at any time t.", [0, 1]],

    # Heterotroph organism paramaters
    "ks": [Decimal('882000000'), "fg C - The organic carbon concentration at which u = 1/2 u0.", [0, 10**10]],
    "organic_carbon_content_per_cell": [Decimal('140'), "fg/cell -  Amount of organic carbon per cell.", [0, 500]],
    "inorganic_carbon_content_per_cell": [Decimal('0'), "fg/cell - Amount of inorganic carbon per cell.", [0, 500]],
    "mu_max": [Decimal('0.06'), "day^-1 - Max growth rate.", [0, 10**5]],
    "growth_penalty": [Decimal('0.5'), "% - Penalty to reduce u_max by. 0 is normal growth.", [0,1]],
    "base_maintenance_per_cell": [Decimal('0.0006') / (Decimal('0.00069') * Decimal('0.4')),
                                  "fg carbon/fg dry mass x day - Constant maintenance energy coefficient.", [0,400]],
    "m_prime": [Decimal('0'),
                "fg glucose/fg dry mass x day - Growth rate dependent maintenance coefficient, known as m_prime in pirt (1982).", [0, 1000]],
    "maximum_growth_yield": [75 / Decimal('72.07'), "fg dry mass/fg carbon x day - Maximum growth yield.", [0, 500]],

    # Inorganic carbon fixing paramaters
    "inorganic_carbon_fixing_factor": [Decimal('0'), "%/day - The percentage of inorganic carbon present fixed.", [0, 1]],

    # System paramaters
    "carrying_capacity": [Decimal('1e9'), "cells/ml - Maximum cell density.", [10**4, 10**10]]
}

default_initial_value = {
    # Initial conditions
    "initial_carbon_content": [Decimal('1.28682e+12'), "fg/ml - Initial amount of organic carbon in the system per ml"],
    "initial_inorganic_carbon_content": [Decimal('0'),
                                         "fg/ml - Initial amount of inorganic carbon in the system per ml."],
    "initial_cell_count": [Decimal('1e5'), "cells/ml - Initial number of cells in the system per ml."],

    # System Paramaters
    "duration": [Decimal('11000') * Decimal('365.25'), "days - Duration of the system."],
}





def model(t, y, paramaters):
    """
    The model.
    """
    ## PARAMATERS
    # Load paramaters
    organic_carbon_input = paramaters["organic_carbon_input"]
    inorganic_carbon_input = paramaters["inorganic_carbon_input"]
    natural_death_fraction = paramaters["natural_death_fraction"]
    ks = paramaters["ks"]
    organic_carbon_content_per_cell = paramaters["organic_carbon_content_per_cell"]
    inorganic_carbon_content_per_cell = paramaters["inorganic_carbon_content_per_cell"]
    mu_max = paramaters["mu_max"]
    growth_penalty = paramaters["growth_penalty"]
    base_maintenance_per_cell = paramaters["base_maintenance_per_cell"]
    m_prime = paramaters["m_prime"]
    maximum_growth_yield = paramaters["maximum_growth_yield"]
    inorganic_carbon_fixing_factor = paramaters["inorganic_carbon_fixing_factor"]
    carrying_capacity = paramaters["carrying_capacity"]

    # Load state conditions
    organic_carbon_content = Decimal(y[0])
    inorganic_carbon_content = Decimal(y[1])
    cell_count = Decimal(y[2])
    t = Decimal(t)


    ## CELL COUNT
    # Growth
    penalized_mu_max = mu_max * (1 - growth_penalty)
    growth = penalized_mu_max * (organic_carbon_content / (ks + organic_carbon_content)) * cell_count * (
            1 - cell_count / carrying_capacity)

    # Specific growth rate
    # next_cell_count = cell_count + growth
    # specific_growth_rate = 0
    # if next_cell_count > 0 and cell_count > 0:
    #     specific_growth_rate = max((np.log(next_cell_count) - np.log(cell_count)), 0)

    # Organic carbon requirement
    required_organic_carbon_per_cell = base_maintenance_per_cell  # + m_prime * (1 - specific_growth_rate/mu_max)
    required_organic_carbon = required_organic_carbon_per_cell * cell_count

    # Starvation
    organic_carbon_missing = max(required_organic_carbon - organic_carbon_content, 0)
    starvation_deaths = 0 if organic_carbon_missing == 0 else organic_carbon_missing / required_organic_carbon_per_cell

    # Natural Death
    natural_deaths = natural_death_fraction * cell_count - natural_death_fraction * starvation_deaths

    # Deaths
    deaths = natural_deaths + starvation_deaths

    # Net cell count change
    dndt = growth - deaths

    ## CARBON
    carbon_consumption = required_organic_carbon_per_cell * (cell_count - deaths)
    fixed_carbon = inorganic_carbon_fixing_factor * inorganic_carbon_content

    # Inorganic carbon
    didt = inorganic_carbon_input + inorganic_carbon_content_per_cell * (
            deaths - growth) + carbon_consumption - fixed_carbon

    # Organic carbon
    dsdt = organic_carbon_input + organic_carbon_content_per_cell * (
            deaths - growth) - carbon_consumption + fixed_carbon

    return [dsdt, didt, dndt]


def plot(S, I, N, t):
    # Organic carbon plot
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(t / 365, S, label='Organic carbon', color="blue")
    plt.xlabel('Years from start')
    plt.ylabel('femtograms/ml')
    plt.title('Organic carbon over time')
    plt.legend(loc=0)

    # Inorganic carbon plot
    plt.subplot(2, 2, 2)
    plt.plot(t / 365, I, label='Inorganic carbon', color="brown")
    plt.xlabel('Years from start')
    plt.ylabel('femtograms/ml')
    plt.title('Inorganic carbon over time')
    plt.legend(loc=0)

    # Cell count plot
    plt.subplot(2, 2, 3)
    plt.plot(t / 365, N, label='Cells', color="green")
    plt.xlabel('Years from start')
    plt.ylabel('cells/ml')
    plt.title('Cell count over time')
    plt.legend(loc=0)

    plt.show()


def calculate_brine_expansion():
    # Lower and upper bounds for how much migration we need. From Go excel file.
    carbon_density_lower = Decimal(51.23) / Decimal('1e6')  # g/ml - Carbon surrounding the brine (lower bound)
    carbon_density_upper = Decimal(1286.81) / Decimal('1e6')  # g/ml - Carbon surrounding the brine (upper bound)

    # Calculate system organic carbon requirements
    total_carbon_required = (I[-1] - initial_carbon_content - organic_carbon_content_per_cell * initial_cell_count) / one_gram_in_femto * brine_volume  # grams
    total_carbon_required_per_year = total_carbon_required / (duration / Decimal('365.25'))  # (grams) /year

    # Calculate volume needed to expand
    volume_needed_lower = total_carbon_required_per_year / carbon_density_lower  # The volume that must be covered to get the required amount of C
    volume_needed_upper = total_carbon_required_per_year / carbon_density_upper  # The volume that must be covered to get the required amount of C

    # Expansion ratio of volime
    ratio_lower = volume_needed_lower / brine_volume
    ratio_upper = volume_needed_upper / brine_volume

    # Expansion ratio of area
    ratio_area_lower = (1 + volume_needed_lower / brine_volume) ** (Decimal('2/3'))
    ratio_area_upper = (1 + volume_needed_upper / brine_volume) ** (Decimal('2/3'))

    # Expansion ratio linearly
    ratio_dimensions_lower = (1 + volume_needed_lower / brine_volume) ** (Decimal('1/3'))
    ratio_dimensions_upper = (1 + volume_needed_upper / brine_volume) ** (Decimal('1/3'))

    a_low = a * ratio_dimensions_lower
    b_low = b * ratio_dimensions_lower

    a_high = a * ratio_dimensions_upper
    b_high = b * ratio_dimensions_upper

    print("Total carbon required over timespan is: " + "{:.3f}".format(total_carbon_required) + "g")
    print("Total carbon required on average per year is: " + "{:.3f}".format(total_carbon_required_per_year) + "g")

    print("Lower bound for % expansion of volume based on brine volume per year is: " + "{:.0f}".format(
        ratio_lower * 100) + "%")
    print("Upper bound for % expansion of volume based on brine volume per year is: " + "{:.0f}".format(
        ratio_upper * 100) + "%")

    print("Lower bound for % expansion of area based on brine volume per year is: " + "{:.0f}".format(
        (ratio_area_lower - 1) * 100) + "%")
    print("Upper bound for % expansion of area based on brine volume per year is: " + "{:.0f}".format(
        (ratio_area_upper - 1) * 100) + "%")

    print("Lower bound of carbon, brine will expand by (a, b) cm per year: (" + "{:.2f}".format(
        a_low) + ", " + "{:.2f}".format(b_low) + ") = " + "{:.0f}".format((ratio_dimensions_lower - 1) * 100) + "%")
    print("Upper bound of carbon, brine will expand by (a, b) cm per year: (" + "{:.2f}".format(
        a_high) + ", " + "{:.2f}".format(b_high) + ") = " + "{:.0f}".format((ratio_dimensions_upper - 1) * 100) + "%")


def sensitivity_analysis_scipy(ivp):
    problem = {
        'num_vars': len(default_paramaters),
        'names': default_paramaters.keys(),
        'bounds': [param[-1] for param in default_paramaters.values()]
    }

    param_values = saltelli.sample(problem, 16384)

    Y = np.zeros([param_values.shape[0]])  # Array to store values in

    for i, paramaters in enumerate(param_values):
        paramaters = {key: Decimal(val) for key, val in zip(default_paramaters.keys(), paramaters)}
        Y[i] = solve_using_scipy(paramaters, ivp).y[0][-1]

    Si = sobol.analyze(problem, Y)
    return Si


def solve_using_scipy(paramaters, ivp):
    y0 = np.asarray([ivp["initial_carbon_content"], ivp["initial_inorganic_carbon_content"], ivp["initial_cell_count"]])
    t_span = (0, ivp["duration"])

    return solve_ivp(fun=model, y0=y0, t_span=t_span, args=[paramaters], method="LSODA", rtol=1e-9, atol=1e-9)


if __name__ == "__main__":
    # Strip help from dicts
    paramaters = {key: val[0] for key, val in zip(default_paramaters.keys(), default_paramaters.values())}
    ivp = {key: val[0] for key, val in zip(default_initial_value.keys(), default_initial_value.values())}

    # Sympy solver
    # sympy_model()

    # Scipy solver
    # soln = solve_using_scipy(paramaters)
    # S, I, N, t = soln.y[0], soln.y[1], soln.y[2], soln.t
    # plot(S, I, N, t)

    # Sensitivity analysis
    Si = sensitivity_analysis_scipy(ivp)
    print(Si)
    axes = Si.plot()
