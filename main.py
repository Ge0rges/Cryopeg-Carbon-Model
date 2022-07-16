'''
This file contains master functions for the different cryopeg brine scenarios.
Run the method that corresponds to the scenario you'd like to simulate.
All relevant paramaters will be calculated for that scenario including:
    - Maintenance energy high and low bounds
    - Model output given the calculated maintenance energy
    - Brine expansion factor
The main function of this file executes all scenarios and the sensitivity analytsis.
Call run_all_analysis() with custom paramaters to create your own scenario.
'''

from maintenance_energy_estimator import *
from brine_expansion_calculator import calculate_brine_expansion
from model import run_model, sensitivity_analysis
from plots import plot_model

def cb1_general_scenario():
    beo_poc = 0
    beo_doc = 0

    # PARAMATERS
    starting_carbon = beo_poc + beo_doc
    ending_carbon = cb1_poc + cb1_doc
    timespan = 40000
    growth_rate = 0.06
    carbon_per_cell =
    inorganic_carbon_per_cell =

    run_all_analysis()

    return

def cb4_scenario():
    return

def cbiw_scenario():
    return

def run_analysis(starting_carbon, ending_carbon, timespan, growth_rate, carbon_per_cell, inorganic_carbon_per_cell):

    # Replace the values in the default dicts

    # Calculate the maintenance energy values
    me_upper =
    me_lower =

    # Replace default maintenance energy value in the dicitonnary

    # Run the model
    run_model()

    # Run brine expansion
    calculate_brine_expansion()

    return

if __name__ == "__main__":
    cb1_general_scenario()
    cb4_scenario()
    cbiw_scenario()
    sensitivy_analysis()
