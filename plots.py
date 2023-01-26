"""
Contains the functions that plot various results using matplotlib and seaborn.
"""

from analysis import Analysis

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

matplotlib.use('TkAgg')
sns.set_theme()


def plot_one_analysis(analysis: Analysis):
    """
    Plots the result of a model iteration: inorganic and organic carbon overlayed, and cell density separately
    in one figure. Returns a figure.
    """
    # Make seaborn data frame
    data = analysis.model_result.get_dataframe(analysis.scenario.title, analysis._variable_title)
    melted = pd.melt(data, id_vars="Years from start", value_vars=("POC", "DOC", "IC", "Cells"), var_name="Variable",
                     value_name="femtograms C/mL or cells/mL")

    # Make figure
    grid = sns.relplot(data=melted, x="Years from start", y="femtograms C/mL or cells/mL", hue="Variable",
                       palette=["brown", "blue", "black", "green"], kind="line")

    grid.set(xscale="log", yscale="log")
    grid.set_titles(template=analysis.title)

    ax = grid.axes.ravel()[0]

    ax.set_ylim([1, 10 ** 14])
    ax.set_xlim([0.01, 10 ** 5])

    grid.tight_layout()

    return grid


def plot_all_scenarios_all_analyses(analyses: [Analysis], color_cycle: int = None):
    """
    Plots the results of many model outputs in one figure. In one figure, creates four columns, one for each data type.
    Each row is an analysis. Returns a figure.
    """

    # Get the right color and dashes for plot decomposition
    cm = sns.color_palette("Paired")
    cp = ["#70381D", cm[7], "#74B6C2"]
    cp = cp if color_cycle is None else [cp[color_cycle]]
    # dashes = [(2, 1), (5, 5), (1, 1), (4, 3)] if len(analyses) <= 4 else None

    # Build the dataframes for each category then melt them
    all_data = pd.DataFrame()
    for analysis in analyses:
        data = analysis.model_result.get_dataframe(analysis.scenario.title, analysis._variable_title)
        all_data = pd.concat([all_data, data]).reset_index(drop=True)

    melted_data = all_data.melt(id_vars=["Years from start", "Scenario", "Analysis type"],
                                value_vars=["POC", "DOC", "IC", "Cells"], var_name="Data type")

    # Plot
    grid = sns.relplot(data=melted_data, x="Years from start", y="value", palette=cp, aspect=0.7,
                       col="Data type", row="Analysis type", hue="Scenario",  # style="Analysis type", dashes=dashes,
                       kind="line", facet_kws={'sharey': False, 'sharex': False, "margin_titles": True,
                                               "legend_out": False})

    grid.set(xscale="log")
    grid.set_titles(template="{col_name} over time", row_template="{row_name}", col_template="{col_name} over time")

    # Tune relplot labels
    y_labels = ["femtograms C/mL"] * 3 + ["cells/mL"]
    y_lim = [[1, 10 ** 14]] * 3 + [[1, 10 ** 10]]
    for i, ax in enumerate(grid.axes.ravel()):
        if i < 4:
            label = y_labels[i % len(y_labels)]
            ax.set_ylabel(label)
            ax.set_xlabel("Years from start")
        else:
            ax.set_ylabel("")
            ax.set_xlabel("")

        lim = y_lim[i % len(y_lim)]
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylim(lim)
        ax.set_xlim([0.01, 10 ** 5])

    grid.tight_layout()
    grid.fig.subplots_adjust(hspace=0.2, wspace=0.4)

    return grid


def plot_multiple_scenarios_one_row(analyses: [Analysis], color_cycle: int = None):
    """
    Plots the results of many model outputs in one figure. In one figure, creates four subplots.
    Each subplot is an overlay of organic carbon, inorganic carbon, and cell densities, from each result.
    Color_cycle shifts the color map to the selected index, used to plot one type of scenario.
    Returns a figure.
    """

    # Get the right color and dashes for plot decomposition
    cm = sns.color_palette("Paired")
    cp = [cm[11], cm[7], cm[1]]
    cp = cp if color_cycle is None else [cp[color_cycle]]
    dashes = [(2, 1), (5, 5), (1, 1), (4, 3)] if len(analyses) <= 4 else None

    # Build the dataframes for each category then melt them
    all_data = pd.DataFrame()
    for analysis in analyses:
        data = analysis.model_result.get_dataframe(analysis.scenario.title, analysis._variable_title)
        all_data = pd.concat([all_data, data]).reset_index(drop=True)

    melted_data = all_data.melt(id_vars=["Years from start", "Scenario", "Analysis type"],
                                value_vars=["POC", "DOC", "IC", "Cells"], var_name="Data type")

    # Plot
    grid = sns.relplot(data=melted_data, x="Years from start", y="value", palette=cp, aspect=0.7,
                       col="Data type", hue="Scenario", style="Analysis type", dashes=dashes,
                       kind="line", facet_kws={'sharey': False, 'sharex': True, "xlim": [0.01, 10 ** 5]})

    grid.set(xscale="log")
    grid.set_titles(template="{col_name} over time")

    # Tune relplot labels
    y_labels = ["femtograms C/mL"] * 3 + ["cells/mL"]
    y_lim = [[1, 10 ** 14]] * 3 + [[1, 10 ** 10]]
    for i, ax in enumerate(grid.axes.ravel()):
        label = y_labels[i % len(y_labels)]
        lim = y_lim[i % len(y_lim)]

        ax.set_ylabel(label)
        ax.set_yscale("log")
        ax.set_ylim(lim)

    return grid.tight_layout()


def plot_one_result_type_all_analyses(analyses: [Analysis], data_type: str, main_index: int = 0, color_cycle: int = None):
    """
    Plots one data type result from each analysis.  Data_type string sets which result type to plot.
    If main_index is set to an index in Analyses, that result is plotted across multiple rows in its own column.
    Color_cycle shifts the color map to the selected index, used to plot one type of scenario.
    Returns a figure.
    """
    # Get the right color and dashes for plot decomposition
    cm = sns.color_palette("Paired")
    cp = ["#70381D", cm[7], "#74B6C2"]
    cp = cp if color_cycle is None else [cp[color_cycle]]

    # Move the main_analysis to front
    analyses.insert(0, analyses.pop(main_index))

    # Build the dataframes for each category then melt them
    all_data = pd.DataFrame()
    for analysis in analyses:
        data = analysis.model_result.get_dataframe(analysis.scenario.title, analysis._variable_title)
        all_data = pd.concat([all_data, data]).reset_index(drop=True)

    melted_data = all_data.melt(id_vars=["Years from start", "Scenario", "Analysis type"],
                                value_vars=[data_type], var_name="Data type")

    # Plot
    grid = sns.relplot(data=melted_data, x="Years from start", y="value", palette=cp, aspect=0.7,
                       col="Analysis type", hue="Scenario", kind="line", col_wrap=4,
                       facet_kws={'sharey': True, 'sharex': True, "legend_out": False})

    grid.set_titles(template="{col_name}")
    grid.set(xscale="log")

    # Tune relplot labels
    y_label = "cells/mL" if data_type == "Cells" else "femtograms C/mL"
    y_lim = [1, 10 ** 10] if data_type == "Cells" else [1, 10 ** 14]
    for i, ax in enumerate(grid.axes.ravel()):
        ax.set_ylabel(y_label)

        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_ylim(y_lim)
        ax.set_xlim([0.01, 10 ** 5])

    grid.tight_layout()
    grid.fig.subplots_adjust(wspace=0.4)

    return grid


def plot_sensitivity(analysis: Analysis):
    """
    Plots the results of a sensitivity analysis as a bar growth per variable. Includes both total and first order
    Sobol indices. Returns a figure.
    """

    # Check there are results
    if analysis.sensitivity_analysis_result is None:
        return None

    # Get dataframes
    df = analysis.sensitivity_analysis_result.get_dataframe(analysis.scenario).set_index('Parameter')

    # We're only plotting mean
    df_mean = df[df.Output.eq('Mean')]

    vals = df_mean[['Total-effect', 'First-order']]
    yerr = df_mean[['Total Error', 'First Error']]
    yerr.columns = vals.columns

    # Use pandas plot function
    ax = vals.plot(kind='bar', yerr=yerr, figsize=(12, 6), ylabel="Value", width=0.8, rot=0, title='Average sensitivity across outputs')

    # Add bar label
    for c in ax.containers:
        if type(c) == matplotlib.container.BarContainer:
            labels = [f'{h:.2f}' if (h := v.get_height()) >= 0.01 else "<0.01" for v in c]
            ax.bar_label(c, labels=labels, label_type='center', padding=5)

    return ax.get_figure()


def hypothetical_growth_scenarios():
    """
    Plots a set of hypothetical growth curves defined by various equations. Returns a figure.
    """

    fig, axis = plt.subplots(1, 1)

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

    fig.tight_layout()

    return fig


if __name__ == "__main__":
    # Generate and save just the hypothetical growth figure
    fig = hypothetical_growth_scenarios()
    fig.savefig("Plots/" + "hypothetical_growth.pdf", format="pdf", dpi=500, bbox_inches='tight')
