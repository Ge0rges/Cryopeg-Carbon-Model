"""
This file contains master functions for the different cryopeg brine scenarios.
Run the method that corresponds to the scenario you'd like to simulate.
All relevant paramaters will be calculated for that scenario including:
    - Maintenance energy high and low bounds
    - Model output given the calculated maintenance energy
    - Brine expansion factor
The main function of this file executes all scenarios and the sensitivity analytsis.
Call run_all_analysis() with custom paramaters to create your own scenario.
"""

import csv

from scenario import *
from plots import *


def log_results(analyses, savelog=True):
    """
    Displays or saves the plots for each analysis done. Prints or save endpoints of model and maintenance energy
    calculations. Saves plots and values to disk if savelog is true.
    """
    csv_header = ["Analysis title", "End dOC", "End cell density", "Maintenance energy lower bounds",
                  "Maintenance energy upper bound", "Minimum growth rate", "Minimum doubling time",
                  "Brine expansion factor"]
    csv_rows = []

    for analysis in analyses:
        # Individual plots & Sensitivity analysis
        model_fig = plot_model(analysis)
        sa_fig = plot_sensitivity(analysis)

        # Endpoints and ME values
        values = [analysis.title,
                  np.format_float_scientific(analysis.model_result.dOC[-1], precision=2),
                  np.format_float_scientific(analysis.model_result.cells[-1], precision=2),
                  np.format_float_scientific(analysis.maintenance_energy_result.lower_bound_me, precision=4),
                  np.format_float_scientific(analysis.maintenance_energy_result.upper_bound_me, precision=4),
                  np.format_float_scientific(analysis.maintenance_energy_result.minimum_growth_rate, precision=4),
                  analysis.maintenance_energy_result.minimum_doubling_time / 365.25,
                  "{:.2f}%".format(analysis.expansion_result.ratio_dimensions)]

        # Save plots and write values to CSV
        if savelog:
            model_fig.savefig("Plots/" + analysis.title + "_model.pdf", format="pdf", dpi=500, bbox_inches='tight')
            csv_rows.append(values)

            if analysis.do_sensitivity_analysis:
                sa_fig.savefig("Plots/" + analysis.title + "_sa.pdf", format="pdf", dpi=500, bbox_inches='tight')

        else:
            model_fig.show()
            if analysis.do_sensitivity_analysis:
                sa_fig.show()

            print(csv_header)
            print(values)

    # Write the values to CSV if required
    if savelog:
        with open('Plots/values.csv', 'w') as f:
            write = csv.writer(f)
            write.writerow(csv_header)
            write.writerows(csv_rows)

    # Make collective plot
    colors = None  #["blue", "orange", "red"]
    line_styles = None #["solid", "dashed", "dotted"]
    scenarios_fig = plot_multiple_scenarios(analyses, colors, line_styles)

    if savelog:  # Plots folder must exist in same directory as main.py
        scenarios_fig.savefig("Plots/all_scenarios_model.pdf", format="pdf", dpi=500, bbox_inches='tight')

    else:  # Show the figures
        scenarios_fig.show()

    return


if __name__ == "__main__":  # Generates all figures and data points.
    # All sensitivity analyses should be the same.
    scenarios = [cb1_scenario(), cb4_scenario(), cbiw_scenario()]
    all_analyses = []

    # On every scenario, try every analysis configuration. Run a sensitivity analysis only on the first one.
    for scenario in scenarios:
        for use_minimum_growth_rate in [True, False]:
            for use_me_lower_bound in [True, False]:
                a = Analysis(scenario, use_minimum_growth_rate, use_me_lower_bound)
                a.run_analysis(use_minimum_growth_rate & use_me_lower_bound)  # Run SA on first config pair only.

                all_analyses.append(a)

    # Save plots and values of all results.
    log_results(all_analyses, savelog=True)

    # Generate the hypothetical growth scenarios figure.
    hg_fig = hypothetical_growth_scenarios()
    hg_fig.savefig("Plots/" + "hypothetical_growth.pdf", format="pdf", dpi=500, bbox_inches='tight')
