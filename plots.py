import matplotlib.pyplot as plt
import numpy as np


def plot_model(S, I, N, t, title):

    fig, axs = plt.subplots(1, 2, constrained_layout=True, dpi=500)

    fig.suptitle(title, fontsize="x-large", fontweight="medium")

    # Carbon plot
    axs[0].loglog(t / 365, S, label='Organic carbon', color="blue", linewidth=2.5)
    axs[0].loglog(t / 365, I, label='Inorganic carbon', color="brown", linewidth=2.5)

    axs[0].set_xlim([0.01, 10**5])
    axs[0].set_xlabel('Years from start')
    axs[0].set_ylabel('femtograms C/mL')
    axs[0].set_title('Carbon over time')
    axs[0].legend(loc=0)

    # Cell count plot
    axs[1].loglog(t / 365, N, label='Cells', color="green", linewidth=2.5)
    axs[1].set_xlabel('Years from start')
    axs[1].set_ylabel('cells/mL')
    axs[1].set_title('Cell count over time')
    axs[1].set_ylim([1, 10**10])
    axs[1].set_xlim([0.01, 10**5])
    axs[1].legend(loc=0)

    return fig


def plot_multiple_scenarios(S_array, I_array, N_array, t_array, labels, title, colors, line_styles):
    fig, axs = plt.subplots(1, 3, dpi=500, figsize=(20, 10))

    fig.suptitle(title, fontsize="xx-large", fontweight="medium")

    # Organic carbon plot
    for label, S, t, color, style in zip(labels, S_array, t_array, colors, line_styles):
        axs[0].loglog(t / 365, S, label=label, color=color, linewidth=2.5, linestyle=style)

    axs[0].set_ylim([0.01, 10**14])
    axs[0].set_xlim([0.01, 10**5])
    axs[0].set_xlabel('Years from start')
    axs[0].set_ylabel('femtograms C/mL')
    axs[0].set_title('Organic carbon over time')

    # Inorganic carbon plot
    for label, I, t, color, style in zip(labels, I_array, t_array, colors, line_styles):
        axs[1].loglog(t / 365, I, label=label, color=color, linewidth=2.5, linestyle=style)

    axs[1].set_ylim([0.01, 10**14])
    axs[1].set_xlim([0.01, 10**5])
    axs[1].set_xlabel('Years from start')
    axs[1].set_ylabel('femtograms C/mL')
    axs[1].set_title('Inorganic carbon over time')

    # Cell count plot
    for label, N, t, color, style in zip(labels, N_array, t_array, colors, line_styles):
        axs[2].loglog(t / 365, N, label=label, color=color, linewidth=2.5, linestyle=style)

    axs[2].set_ylim([1, 10**10])
    axs[2].set_xlim([0.01, 10**5])
    axs[2].set_xlabel('Years from start')
    axs[2].set_ylabel('cells/mL')
    axs[2].set_title('Cell count over time')

    axs[2].legend(loc=0)

    return fig


def hypothetical_growth_scenarios():
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


def plot_sensitivity(ST, S1, p_names, label):
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

    axs.set_xlabel('Organism parameter', fontsize=18)
    axs.set_ylabel("Sobol Index", fontsize=18)
    axs.set_xticks([r + barWidth for r in range(len(ST))], p_names)
    axs.set_title("Sensitivity analysis of organism parameters - " + label + " scenario")

    axs.legend(loc=0)

    return fig


if __name__ == "__main__":
    fig = hypothetical_growth_scenarios()
    fig.savefig("Plots/" + "hypothetical_growth.pdf", format="pdf", dpi=500, bbox_inches='tight')
