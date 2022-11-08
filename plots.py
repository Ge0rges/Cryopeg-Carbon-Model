"""
Contains the functions that plot various results using matplotlib and seaborn.
"""

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from analysis import Analysis
import pandas as pd

matplotlib.use('TkAgg')
sns.set_theme()


def plot_model(analysis: Analysis):
    """
    Plots the result of a model iteration: inorganic and organic carbon overlayed, and cell density seperately
    in one figure. Returns a figure.
    """
    # Make seaborn data frame
    data = analysis.model_result.get_dataframe(analysis.scenario.title, analysis._variable_title)
    melted = pd.melt(data, id_vars=("Years from start"), value_vars=("POC", "DOC", "IC"), var_name=("Carbon type"),
                     value_name=("femtograms C/mL"))

    # Make figure
    f, axs = plt.subplots(1, 2, figsize=(8, 4))
    sns.lineplot(data=melted, x="Years from start", y="femtograms C/mL", hue="Carbon type",
                 palette=['brown', 'blue', 'black'], ax=axs[0])
    sns.lineplot(data=data, x="Years from start", y="Cells", color="green", ax=axs[1])

    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')

    axs[0].set_ylim([1, 10 ** 14])
    axs[0].set_xlim([0.01, 10 ** 5])
    axs[1].set_ylim([1, 10 ** 10])
    axs[1].set_xlim([0.01, 10 ** 5])

    axs[1].set_ylabel("cells/mL")

    axs[0].set_title('Carbon over time')
    axs[1].set_title('Cells over time')

    f.tight_layout()

    return f


def plot_multiple_scenarios(analyses: [Analysis], advance_cycler: int = None):
    """
    Plots the results of many model outputs in one figure. In one figure, creates three subplots.
    Each subplot is an overlay of organic carbon, inorganic carbon, and cell densities, from each result.
    Returns a figure.
    """

    # Get the right color for plot decomposition
    cm = sns.color_palette("Paired")
    cp = [cm[11], cm[7], cm[1]]
    cp = cp if advance_cycler is None else [cp[advance_cycler]]

    # Build the dataframes for each category then melt them
    all_data = pd.DataFrame()
    for analysis in analyses:
        data = analysis.model_result.get_dataframe(analysis.scenario.title, analysis._variable_title)
        all_data = pd.concat([all_data, data]).reset_index(drop=True)

    melted_data = all_data.melt(id_vars=["Years from start", "Scenario", "Analysis type"],
                                value_vars=["POC", "DOC", "IC", "Cells"], var_name=("Data type"))

    # Plot
    grid = sns.relplot(data=melted_data, x="Years from start", y="value", palette=cp, aspect=0.7,
                       col="Data type", hue="Scenario", style="Analysis type", dashes=[(2, 1), (5, 5)],
                       kind="line", facet_kws={'sharey': False, 'sharex': True, "xlim": [0.01, 10 ** 5]})

    grid.set(xscale="log", yscale="log")
    grid.set_titles(template="{col_name} over time")

    # Tune relplot
    y_labels = ["femtograms C/mL"] * 3 + ["cells/mL"]
    y_lim = [[1, 10 ** 14]] * 3 + [[1, 10 ** 10]]
    for ax, label, lim in zip(grid.axes.ravel(), y_labels, y_lim):
        ax.set_ylabel(label)
        ax.set_yscale("log")
        ax.set_ylim(lim)

    return grid.tight_layout()


def plot_sensitivity(analysis: Analysis):
    """
    Plots the results of a sensitivity analysis as a bar growth per variable. Includes both total and first order
    Sobol indices. Returns a figure.
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
    barWidth = 0.28
    br1 = np.arange(len(ST))
    br2 = [x + barWidth for x in br1]

    # Total sobol index
    axs.bar(br1, ST, width=barWidth, label="Total-effect")
    axs.bar(br2, S1, width=barWidth, label="First-order")

    for index in range(len(br1)):
        axs.text(br1[index] - barWidth/2, ST[index]+0.01, "%.2f" % ST[index], size=12)
        axs.text(br2[index] - barWidth/2, S1[index]+0.01, "%.2f" % S1[index], size=12)

    axs.set_xlabel('Organism parameter', fontsize=12)
    axs.set_ylabel("Sobol Index", fontsize=12)
    axs.set_xticks([r + barWidth for r in range(len(ST))], p_names)

    axs.legend(loc=0)

    return fig


def hypothetical_growth_scenarios():
    """
    Plots a set of hypothetical growth curves defined by various equations. Returns a figure.
    """

    fig, axis = plt.subplots(1, 1, tight_layout=True, dpi=500)

    x = np.linspace(0, 40000, num=4000000)

    x_exp = np.linspace(0, 6.90776, num=len(x))
    y_exp = 10**5 * np.exp(x_exp)
    y_cyclic = 10**5 + (10**8 - 10**5)/2 + np.sin(x/4244.1333333333 + 1.5*np.pi) * (10**8 - 10**5)/2
    y_ng = np.full(shape=len(x), fill_value=10**8)
    y_spike = np.where(x < 0.02, 10**5, np.where(x < 10, x*10**7 + 10*5, 10**8))

    axis.loglog(x, y_spike, label='Rapid', color="blue")
    axis.loglog(x, y_exp, label='Slow', color="green")
    axis.loglog(x, y_cyclic, label='Carbon addition', color="red")
    axis.loglog(x, y_ng, label='No Growth', color="brown", linestyle='dashed')

    axis.set_ylim([0.1, 10**10])
    axis.set_xlim([0.01, 10**5])
    axis.set_xlabel('Years from start')
    axis.set_ylabel('cells/mL')
    axis.legend(loc=4)

    return fig


if __name__ == "__main__":
    # Generate and save just the hypothetical growth figure
    fig = hypothetical_growth_scenarios()
    fig.savefig("Plots/" + "hypothetical_growth.pdf", format="pdf", dpi=500, bbox_inches='tight')
