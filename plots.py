import matplotlib.pyplot as plt


def plot_model(S, I, N, t):
    # Organic carbon plot
    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(t / 365, S, label='Organic carbon', color="blue")
    plt.xlabel('Years from start')
    plt.ylabel('femtograms/ml')
    plt.title('Organic carbon over time')
    plt.legend(loc=0)

    # Inorganic carbon plot
    plt.subplot(2, 2, 2)
    plt.plot(t / 365, I, label='Inorganic carbon', color="brown")
    plt.xlabel('Years from start')
    plt.ylabel('femtograms/ml')
    plt.title('Inorganic carbon over time')
    plt.legend(loc=0)

    # Cell count plot
    plt.subplot(2, 2, 3)
    plt.plot(t / 365, N, label='Cells', color="green")
    plt.xlabel('Years from start')
    plt.ylabel('cells/ml')
    plt.title('Cell count over time')
    plt.legend(loc=0)

    plt.show()


def plot_sensitivity(ST, S1, p_names):
    return
