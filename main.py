"""
This file contains master functions for the different cryopeg brine scenarios.
Run the method that corresponds to the scenario you'd like to simulate.
All relevant parameters will be calculated for that scenario including:
    - Maintenance energy high and low bounds
    - Model output given the calculated maintenance energy
    - Brine expansion factor
The main function of this file executes all scenarios and the sensitivity analysis.
Call run_all_analysis() with custom parameters to create your own scenario.
"""

import csv
import os

from scenario import *
from plots import *


def log_results(analyses):
    """
    Displays or saves the plots for each analysis done. Prints or save endpoints of model and maintenance energy
    calculations. Saves plots and values to disk if savelog is true.
    """
    csv_header = ["Analysis title", "Surrounding pOC", "Brine pOC", "Predicted brine pOC", "Surrounding dOC",
                  "Brine dOC", "Predicted brine dOC", "Present cell density", "Predicted end cell density",
                  "EEA Lower", "EEA Upper", "Real EEA", "EEA predicted timespan", "dOC/cell",
                  "Maintenance energy lower bounds", "Maintenance energy upper bound", "Simulation growth rate",
                  "Minimum growth rate", "Minimum doubling time", "Growth yield", "Brine expansion factor"]
    csv_rows = []

    for analysis in analyses:
        # Individual plots & Sensitivity analysis
        model_fig = plot_one_analysis(analysis)
        sa_fig = plot_sensitivity(analysis)

        # Endpoints and ME values
        values = [analysis.title,
                  # POC
                  np.format_float_scientific(analysis.scenario.start_poc, precision=2),
                  np.format_float_scientific(analysis.scenario.end_poc, precision=2),
                  np.format_float_scientific(analysis.model_result.pOC[-1], precision=2),

                  # DOC
                  np.format_float_scientific(analysis.scenario.start_doc, precision=2),
                  np.format_float_scientific(analysis.scenario.end_doc, precision=2),
                  np.format_float_scientific(analysis.model_result.dOC[-1], precision=2),

                  # Cell density
                  np.format_float_scientific(analysis.scenario.observed_end_cell_density, precision=2),
                  np.format_float_scientific(analysis.model_result.cells[-1], precision=2),

                  # EEA
                  np.format_float_scientific(analysis.eea_estimation.eea_lower, precision=2),
                  np.format_float_scientific(analysis.eea_estimation.eea_upper, precision=2),
                  np.format_float_scientific(analysis.scenario._eea_rate, precision=2),
                  np.format_float_scientific(analysis.eea_estimation.predicted_timespan/365.25, precision=2),

                  # Cell carbon content
                  np.format_float_scientific(analysis.scenario.dissolved_organic_carbon_per_cell, precision=2),

                  # Maintenance energy
                  np.format_float_scientific(analysis.maintenance_energy_result.lower_bound_me, precision=4),
                  np.format_float_scientific(analysis.maintenance_energy_result.upper_bound_me, precision=4),

                  # Growth rate
                  np.format_float_scientific(analysis.scenario._growth_rate, precision=4),
                  np.format_float_scientific(analysis.maintenance_energy_result.minimum_growth_rate, precision=4),
                  np.format_float_scientific(analysis.maintenance_energy_result.minimum_doubling_time / 365.25, precision=2),

                  # Growth yield
                  np.format_float_scientific(analysis.growth_yield, precision=4),

                  # Brine expansion
                  np.format_float_scientific(analysis.expansion_result.ratio_dimensions, precision=2)]

        # Save plots and write values to CSV
        # Make a plots folder is it doesn't exist
        if not os.path.exists("Results/"):
            os.mkdir('Results/')

        model_fig.savefig("Results/" + analysis.title + "_model.pdf", format="pdf", dpi=500)
        csv_rows.append(values)

        if sa_fig:
            sa_fig.savefig("Results/" + analysis.title + "_sa.pdf", format="pdf", dpi=500)

    # Write the values to CSV if required
    with open('Results/values.csv', 'w+') as f:
        write = csv.writer(f)
        write.writerow(csv_header)
        write.writerows(csv_rows)

    # # Plot for selected analyses
    # selected_analyses = []
    # for analysis in analyses:
    #     if any(map(analysis.title.__contains__, ["Min μ - High ME - Measured EEA", "Max μ - Low ME - Measured EEA"])):
    #         selected_analyses.append(analysis)
    #
    # sc_index = int(len(selected_analyses) / len(scenarios))
    # individual_scenarios = [plot_multiple_scenarios_one_row(selected_analyses[(i-1)*sc_index:sc_index*i], color_cycle=i-1) for i in range(1, len(scenarios)+1)]

    # Plot for all scenarios all analyses
    all_analyses_fig = plot_all_scenarios_all_analyses(analyses)
    single_datatype_plot_cells = plot_one_result_type_all_analyses(analyses, "Cells", 0)

    # for i, fig in enumerate(individual_scenarios):
    #     fig.savefig("Results/pair"+str(i)+"_scenarios_model.pdf", format="pdf", dpi=500)

    all_analyses_fig.savefig("Results/all_model_outputs.pdf", format="pdf", dpi=500)
    single_datatype_plot_cells.savefig("Results/all_model_outputs_just_cells.pdf", format="pdf", dpi=500)

    return


if __name__ == "__main__":  # Generates all figures and data points.
    # All sensitivity analyses should be the same.
    scenarios = [cb1_scenario(), cb4_scenario(), cbiw_scenario()]
    all_analyses = []

    # On every scenario, try every analysis configuration. Run a sensitivity analysis only on the first one.
    do_SA = True
    for use_minimum_growth_rate in [True, False]:
        for use_me_lower_bound in [True, False]:
            for use_eea_average in [True, False]:
                for scenario in scenarios:
                    a = Analysis(scenario, use_minimum_growth_rate, use_me_lower_bound, use_eea_average)

                    # Run SA on first config pair - first scenario, only.
                    a.run_analysis(do_sensitivity_analysis=do_SA)
                    do_SA = False

                    all_analyses.append(a)

    # Save plots and values of all results.
    log_results(all_analyses)

    # Generate the hypothetical growth scenarios figure.
    hg_fig = hypothetical_growth_scenarios()
    hg_fig.savefig("Results/" + "hypothetical_growth.pdf", format="pdf", dpi=500)
