from plots import plot_model
from julia import Main
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
default_paramaters = [10 ** 9, 0, 0, 0.001, 0, 882000000, 140, 0, 0.06, 0.0006/(0.00069 * 0.4), "empty"]
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

if __name__ == "__main__":
    S, I, N, t = Main.run_model(default_paramaters, default_ivp)
    plot_model(S, I, N, t)
    # sa_res = Main.run_sensitivity_analysis(default_paramaters, default_paramater_bounds, default_ivp)
    # print(sa_res)
