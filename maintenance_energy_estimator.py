import sympy as sp
import numpy as np
from sympy import symbols, exp, integrate

start_cell = 10 ** 5  # cells/ml - Start cell abundance per ml
end_cell = 10 ** 8  # cells/ml - End cell abundance per ml
timespan = 11000 * 365.25  # days - Community lifespan
start_carbon = 100000  # fg C/ml - Start carbon concentration
end_carbon = 50000  # fg C/ml - End carbon concentration
cell_carbon_content = 400  # fg C - Carbon content a cell takes up/releases upon birth/death


def estimate_no_growth():
    m = (end_carbon - start_carbon)/(end_cell*timespan)

    print("For growth rate 0, maintenance energy value is: " + str(m) + " fg C/(cell * day)")


def estimate_exponential_growth():
    sp.init_printing()

    # N(t) = N_0 * exp(μ*t) where μ is growth rate
    N0, mu, t = symbols("N_0 mu t", real=True)

    N0 = start_cell
    mu = np.log(end_cell / start_cell) / timespan  # days ^-1
    # number_of_generations = np.log(end_cell / start_cell) / np.log(2)  # generations
    # doubling_time = timespan/number_of_generations  # days

    N = N0 * exp(mu*t)

    # dC/dt = -g*dN/dt - Nm   where: g is cell_carbon_content,
    # m is maintenance energy, C is carbon concentration, N is cell concentration
    # Integrate from 0 to timespan.
    # Cf - C0 = -g(Nf-N0) - m∫N(t)dt
    # m = (g(Nf-N0) + Cf - C0)/∫N(t)dt
    m = (cell_carbon_content * (end_cell - start_cell) + end_carbon - start_carbon)/integrate(N, (t, 0, timespan))

    print("For strictly exponential growth, growth rate is: " + str(mu) + "/days")
    print("Maintenance energy value is: " + str(m) + " fg C/(cell * day)")


if __name__ == "__main__":
    estimate_no_growth()
    estimate_exponential_growth()
