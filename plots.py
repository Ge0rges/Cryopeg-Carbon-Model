"""
Contains the functions that plot various results using matplotlib.
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from analysis import Analysis

matplotlib.use('TkAgg')


def plot_model(analysis: Analysis):
    """
    Plots the result of a model iteration: inorganic and organic carbon overlayed, and cell density seperately
    in one figure. Gives the figure a title and returns it.
    """

    fig, axs = plt.subplots(1, 2, constrained_layout=True, dpi=500)

    fig.suptitle(analysis.title, fontsize="x-large", fontweight="medium")

    t = analysis.model_result.t

    # Carbon plot
    axs[0].loglog(t / 365, analysis.model_result.pOC, label='Particulate organic carbon', color="brown", linewidth=2.5)
    axs[0].loglog(t / 365, analysis.model_result.dOC, label='Dissolved organic carbon', color="blue", linewidth=2.5)
    axs[0].loglog(t / 365, analysis.model_result.IC, label='Inorganic carbon', color="black", linewidth=2.5)

    axs[0].set_xlim([0.01, 10**5])
    axs[0].set_xlabel('Years from start')
    axs[0].set_ylabel('femtograms C/mL')
    axs[0].set_title('Carbon over time')
    axs[0].legend(loc=0)

    # Cell count plot
    axs[1].loglog(t / 365, analysis.model_result.cells, label='Cells', color="green", linewidth=2.5)
    axs[1].set_xlabel('Years from start')
    axs[1].set_ylabel('cells/mL')
    axs[1].set_title('Cell count over time')
    axs[1].set_ylim([1, 10**10])
    axs[1].set_xlim([0.01, 10**5])
    axs[1].legend(loc=0)

    return fig


def plot_multiple_scenarios(analyses: [Analysis], colors, line_styles):
    """
    Plots the results of many model outputs in one figure. In one figure, creates three subplots.
    Each subplot is an overlay of organic carbon, inorganic carbon, and cell densities, from each result.
    labels the figure, customizable color and line style for each result set. Returns the figure.
    """
    P_array = [analysis.model_result.pOC for analysis in analyses]
    D_array = [analysis.model_result.dOC for analysis in analyses]
    I_array = [analysis.model_result.IC for analysis in analyses]
    N_array = [analysis.model_result.cells for analysis in analyses]
    t_array = [analysis.model_result.t for analysis in analyses]
    labels = [analysis.title for analysis in analyses]

    assert len(labels) == len(P_array) == len(t_array) == len(D_array) == len(I_array) == len(N_array)

    fig, axs = plt.subplots(1, 4, dpi=500, figsize=(20, 10))

    fig.suptitle("Model outputs of all considered scenarios", fontsize="xx-large", fontweight="medium")

    # Particulate organic carbon plot
    for label, P, t in zip(labels, P_array, t_array):  #, colors, line_styles):
        axs[0].loglog(t / 365, P, label=label, linewidth=2.5)  # color=color, linestyle=style)

    axs[0].set_ylim([0.01, 10**14])
    axs[0].set_xlim([0.01, 10**5])
    axs[0].set_xlabel('Years from start')
    axs[0].set_ylabel('femtograms pOC/mL')
    axs[0].set_title('Particulate organic carbon over time')

    # Dissolved organic carbon plot
    for label, D, t in zip(labels, D_array, t_array):  #, colors, line_styles):
        axs[1].loglog(t / 365, D, label=label, linewidth=2.5)  #, color=color,linestyle=style)

    axs[1].set_ylim([0.01, 10**14])
    axs[1].set_xlim([0.01, 10**5])
    axs[1].set_xlabel('Years from start')
    axs[1].set_ylabel('femtograms dOC/mL')
    axs[1].set_title('Dissolved organic carbon over time')

    # Inorganic carbon plot
    for label, I, t in zip(labels, I_array, t_array):  #, colors, line_styles):
        axs[2].loglog(t / 365, I, label=label, linewidth=2.5)  #, color=color, linestyle=style)

    axs[2].set_ylim([0.01, 10**14])
    axs[2].set_xlim([0.01, 10**5])
    axs[2].set_xlabel('Years from start')
    axs[2].set_ylabel('femtograms C/mL')
    axs[2].set_title('Inorganic carbon over time')

    # Cell count plot
    for label, N, t in zip(labels, N_array, t_array):  #, colors, line_styles):
        axs[3].loglog(t / 365, N, label=label, linewidth=2.5)  #, color=color,linestyle=style)

    axs[3].set_ylim([1, 10**10])
    axs[3].set_xlim([0.01, 10**5])
    axs[3].set_xlabel('Years from start')
    axs[3].set_ylabel('cells/mL')
    axs[3].set_title('Cell count over time')

    axs[2].legend(loc=0)

    return fig


def hypothetical_growth_scenarios():
    """
    Plots a set of hypothetical growth curves defined by various equations. Returns a figure.
    """

    fig, axis = plt.subplots(1, 1, tight_layout=True, dpi=500)

    fig.suptitle("Hypothetical growth scenarios", fontsize="x-large", fontweight="medium")

    x = np.linspace(0, 40000, num=4000000)

    x_exp = np.linspace(0, 6.90776, num=len(x))
    y_exp = 10**5 * np.exp(x_exp)
    y_cyclic = 10**5 + (10**8 - 10**5)/2 + np.sin(x/4244.1333333333 + 1.5*np.pi) * (10**8 - 10**5)/2
    y_ng = np.full(shape=len(x), fill_value=10**8)
    y_spike = np.where(x < 0.02, 10**5, np.where(x < 10, x*10**7 + 10*5, 10**8))

    axis.loglog(x, y_spike, label='Rapid', color="blue")
    axis.loglog(x, y_exp, label='Slow', color="green")
    axis.loglog(x, y_cyclic, label='Carbon addition', color="red")
    axis.loglog(x, y_ng, label='No Growth', color="yellow", linestyle='dashed')

    axis.set_ylim([0.1, 10**10])
    axis.set_xlim([0.01, 10**5])
    axis.set_xlabel('Years from start')
    axis.set_ylabel('cells/mL')
    axis.legend(loc=4)

    return fig


def plot_sensitivity(analysis: Analysis):
    """
    Plots the results of a sensitivity analysis as a bar growth per variable. Includes both total and first order
    Sobol indices. Titles the figure and returns it.
    """

    # Check there are results
    if analysis.sensitivity_analysis_result is None:
        return None

    p_names = analysis.scenario._paramater_names
    ST = analysis.sensitivity_analysis_result.total_sobol_indices
    S1 = analysis.sensitivity_analysis_result.first_order_sobol_indices

    # Filter out zeros (analysis wasn't run)
    indices = []
    for i, x in enumerate(p_names):
        if ST[i] != 0:
            indices.append(i)

    ST = ST[indices]
    S1 = S1[indices]
    p_names = np.asarray(p_names)[indices]

    # Make plot
    fig, axs = plt.subplots(1, 1, tight_layout=True, dpi=500)

    # Set position of bar on X axis
    barWidth = 0.25
    br1 = np.arange(len(ST))
    br2 = [x + barWidth for x in br1]

    # Total sobol index
    axs.bar(br1, ST, width=barWidth, label="Total-effect")
    axs.bar(br2, S1, width=barWidth, label="First-order")

    for index in range(len(br1)):
        axs.text(br1[index] - barWidth/2, ST[index], "%.2f" % ST[index], size=12)
        axs.text(br2[index] - barWidth/2, S1[index], "%.2f" % S1[index], size=12)

    axs.set_xlabel('Organism parameter', fontsize=12)
    axs.set_ylabel("Sobol Index", fontsize=12)
    axs.set_xticks([r + barWidth for r in range(len(ST))], p_names)
    axs.set_title("Sensitivity analysis of organism parameters - " + analysis.scenario.title + " scenario")

    axs.legend(loc=0)

    return fig


if __name__ == "__main__":
    # Generate and save just the hypothetical growth figure
    fig = hypothetical_growth_scenarios()
    fig.savefig("Plots/" + "hypothetical_growth.pdf", format="pdf", dpi=500, bbox_inches='tight')
