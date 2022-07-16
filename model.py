import numpy as np
from diffeqpy import de
from SALib.sample import saltelli
from SALib.analyze import sobol
from plots import plot_model
import matplotlib.pyplot as plt

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
# # Organic carbon around brine
# utgiagvik_nutrients_per_cm2_for_100cm = Decimal('83.4') / one_m2_in_cm2 * one_kilo_in_femto  # fg/cm2 - Initial C available per cm2 over 100cm depth if SOCC100.
# utgiagvik_nutrients_depth = 100  # cm - The NCSCD gives cm2 numbers over a total depth.

default_paramaters = [10 ** 9, 0, 0, 0.001, 0, 882000000, 140, 0, 0.06, 0.0006/(0.00069 * 0.4), None]
default_paramater_names = ["carrying_capacity", "organic_carbon_input", "inorganic_carbon_input",
                           "Death %", "IC F%", "Ks",
                           "OC/Cell", "IC/Cell", "mu_max",
                           "ME", "m_prime"]
default_paramater_bounds = [None, None, None, None, [0, 1], [0, 10**10], [0, 500], [0, 500], None, [0, 400], None]
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


default_ivp = [1286820000000, 0, 10 ** 5, 11000 * 365.25]
default_ivp_names = ["initial_carbon_content", "initial_inorganic_carbon_content", "initial_cell_count", "duration"]

# initial_carbon_content in fg/ml - Initial amount of organic carbon in the system per ml
# initial_inorganic_carbon_content in fg/ml - Initial amount of inorganic carbon in the system per ml.
# initial_cell_count in cells/ml - Initial number of cells in the system per ml.
# duration in days - Duration of the system.



def model(du, u, paramaters, t):
    """
    The model.
    """
    # PARAMATERS
    # Load paramaters - Can't pass dict to julia
    carrying_capacity = paramaters[0]
    organic_carbon_input = paramaters[1]
    inorganic_carbon_input = paramaters[2]
    natural_death_fraction = paramaters[3]
    inorganic_carbon_fixing_factor = paramaters[4]
    ks = paramaters[5]
    organic_carbon_content_per_cell = paramaters[6]
    inorganic_carbon_content_per_cell = paramaters[7]
    mu_max = paramaters[8]
    base_maintenance_per_cell = paramaters[9]
    m_prime = paramaters[10]

    assert m_prime is None  # m_prime not currently taken into account.

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
    # required_organic_carbon_per_cell = base_maintenance_per_cell  + m_prime * (1 - specific_growth_rate/mu_max)

    # Organic carbon requirement
    required_organic_carbon_per_cell = base_maintenance_per_cell
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


def run_model(paramaters, ivp):
    # Create ODE problem
    y0 = ivp[:-1]
    t_span = np.asarray([0, ivp[-1]])
    prob = de.ODEProblem(model, y0, t_span, paramaters)

    # Julia solver
    solution = de.solve(prob, reltol=1e-9, abstol=1e-9)

    S, I, N, t = [x[0] for x in solution.u], [x[1] for x in solution.u], [x[2] for x in solution.u], solution.t
    return S, I, N, t


def sensitivity_analysis(paramaters, paramater_bounds, paramater_names, ivp):
    # Filter out paramaters with no bounds
    p_bounds = []
    p_names = []
    add_back = []

    for i, x in enumerate(paramater_bounds):
        add_back.append("fill")
        if x is None:
            add_back[i] = paramaters[i]
        else:
            p_bounds.append(paramater_bounds[i])
            p_names.append(paramater_names[i])

    # Define the problem inputs for SALib
    problem = {'num_vars': len(p_bounds), 'names': p_names, 'bounds': p_bounds}

    # Sample paramater space for Sobol
    param_values = saltelli.sample(problem, 32)

    # Run the model on each paramater set
    results = np.zeros([param_values.shape[0]])

    for i, params in enumerate(param_values):
        # Fill in Nones with paramaters
        fixed_parmater = [x for x in add_back]
        j = 0
        for k, x in enumerate(fixed_parmater):
            if x == "fill":
                fixed_parmater[k] = params[j]
                j += 1

        S, I, N, t = run_model(fixed_parmater, ivp)
        results[i] = N[-1]

    # Sobol analysis
    Si = sobol.analyze(problem, results)

    print(Si['S1'])
    print(Si['ST'])
    print(p_names)

    ax1, ax2, ax3 = Si.plot()
    plt.show()

    return Si


if __name__ == "__main__":
    # S, I, N, t = run_model(default_paramaters, default_ivp)
    # plot_model(S, I, N, t)
    sensitivity_analysis(default_paramaters, default_paramater_bounds, default_paramater_names, default_ivp)

# def sensitivity_analysis(paramaters, ivp):
#     # Get an array of paramater bounds
#     p_range = [x[-1] if x[-1] is not None else [x[0], x[0]] for x in paramaters.values()]
#
#     # Strip help string and bounds from dicts
#     paramaters = {key: val[0] for key, val in zip(paramaters.keys(), paramaters.values())}
#     ivp = {key: val[0] for key, val in zip(ivp.keys(), ivp.values())}
#
#     # Create ODE problem
#     y0 = np.asarray([ivp["initial_carbon_content"], ivp["initial_inorganic_carbon_content"], ivp["initial_cell_count"]])
#     t_span = np.asarray([0, ivp["duration"]])
#
#     # This function remakes the problem and solves it using new paramater set.
#     def f1(p):
#         prob = de.ODEProblem(model, y0, t_span, p.values())
#         sol = de.solve(prob, reltol=1e-9, abstol=1e-9)
#         return [sol[1, :][-1]]
#
#     de.GlobalSensitivity.gsa(f1, de.GlobalSensitivity.Sobol(), p_range, N=1000)
#
#     return
