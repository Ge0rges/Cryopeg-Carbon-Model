import matplotlib.pyplot as plt
import numpy as np
from diffeqpy import de

from decimal import Decimal, getcontext, ROUND_CEILING

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

## MODEL PARAMATERS
# # Organic carbon around brine
# utgiagvik_nutrients_per_cm2_for_100cm = Decimal('83.4') / one_m2_in_cm2 * one_kilo_in_femto  # fg/cm2 - Initial C available per cm2 over 100cm depth if SOCC100.
# utgiagvik_nutrients_depth = 100  # cm - The NCSCD gives cm2 numbers over a total depth.

default_paramaters = {
    # Inputs
    "organic_carbon_input": [Decimal('0'), "fg/(ml x day) - Organic arbon input per day.", [0, 10 ** 10]],
    "inorganic_carbon_input": [Decimal('0'), "fg/(ml x day) - Inorganic carbon input per day.", [0, 10 ** 10]],

    # Community paramaters
    "natural_death_fraction": [Decimal('0.001'), "%/day - Fraction of cells present that will die at any time t.",
                               [0, 1]],

    # Heterotroph organism paramaters
    "ks": [Decimal('882000000'), "fg C - The organic carbon concentration at which u = 1/2 u0.", [0, 10 ** 10]],
    "organic_carbon_content_per_cell": [Decimal('140'), "fg/cell -  Amount of organic carbon per cell.", [0, 500]],
    "inorganic_carbon_content_per_cell": [Decimal('0'), "fg/cell - Amount of inorganic carbon per cell.", [0, 500]],
    "mu_max": [Decimal('0.06'), "day^-1 - Max growth rate.", [0, 10 ** 5]],
    "base_maintenance_per_cell": [Decimal('0.0006') / (Decimal('0.00069') * Decimal('0.4')),
                                  "fg carbon/fg dry mass x day - Constant maintenance energy coefficient.", [0, 400]],
    "m_prime": [Decimal('0'),
                "fg glucose/fg dry mass x day - Growth rate dependent maintenance coefficient, known as m_prime in pirt (1982).",
                [0, 1000]],
    "maximum_growth_yield": [75 / Decimal('72.07'), "fg dry mass/fg carbon x day - Maximum growth yield.", [0, 500]],

    # Inorganic carbon fixing paramaters
    "inorganic_carbon_fixing_factor": [Decimal('0'), "%/day - The percentage of inorganic carbon present fixed.",
                                       [0, 1]],

    # System paramaters
    "carrying_capacity": [Decimal('1e9'), "cells/ml - Maximum cell density.", [10 ** 4, 10 ** 10]]
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


def model(du, u, paramaters, t):
    """
    The model.
    """
    # PARAMATERS
    # Load paramaters
    organic_carbon_input = paramaters[0]
    inorganic_carbon_input = paramaters[1]
    natural_death_fraction = paramaters[2]
    ks = paramaters[3]
    organic_carbon_content_per_cell = paramaters[4]
    inorganic_carbon_content_per_cell = paramaters[5]
    mu_max = paramaters[6]
    base_maintenance_per_cell = paramaters[7]
    m_prime = paramaters[8]
    maximum_growth_yield = paramaters[9]
    inorganic_carbon_fixing_factor = paramaters[10]
    carrying_capacity = paramaters[11]

    # Load state conditions
    organic_carbon_content = u[0]
    inorganic_carbon_content = u[1]
    cell_count = u[2]


    ## CELL COUNT
    # Growth
    growth = mu_max * (organic_carbon_content / (ks + organic_carbon_content)) * cell_count * (
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
    du[2] = growth - deaths

    ## CARBON
    carbon_consumption = required_organic_carbon_per_cell * (cell_count - deaths)
    fixed_carbon = inorganic_carbon_fixing_factor * inorganic_carbon_content

    # Inorganic carbon
    du[1] = inorganic_carbon_input + inorganic_carbon_content_per_cell * (
            deaths - growth) + carbon_consumption - fixed_carbon

    # Organic carbon
    du[0] = organic_carbon_input + organic_carbon_content_per_cell * (
            deaths - growth) - carbon_consumption + fixed_carbon




def solve_using_julia(paramaters, ivp):
    y0 = np.asarray([float(ivp["initial_carbon_content"]), float(ivp["initial_inorganic_carbon_content"]),
                     float(ivp["initial_cell_count"])])
    t_span = np.asarray([0, float(ivp["duration"])])
    prob = de.ODEProblem(model, y0, t_span, np.asarray([float(x) for x in paramaters.values()]))

    return de.solve(prob, reltol=1e-9, abstol=1e-9)


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


if __name__ == "__main__":
    # Strip help string and bounds from dicts
    paramaters = {key: val[0] for key, val in zip(default_paramaters.keys(), default_paramaters.values())}
    ivp = {key: val[0] for key, val in zip(default_initial_value.keys(), default_initial_value.values())}

    # Julia solver
    soln = solve_using_julia(paramaters, ivp)
    S, I, N, t = [x[0] for x in soln.u],  [x[1] for x in soln.u], [x[2] for x in soln.u], soln.t

    plot(S, I, N, t)

    # # Sensitivity analysis
