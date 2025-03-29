import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter
import json

def plot_inhomo_log(N, Alpha_values, filename):

    with open(filename, 'r') as f:

            data = json.load(f)

    marker = [r'$\otimes$', r'$\ast$', '+', 'x', r'$\star$', r'$\circ$']
    markersize = [3, 4, 4, 3.5, 4, 3]
    markeredgewidth = [0.25, 0.25, 0.3, 0.3, 0.2, 0.05]
    labels = [r'uniform', r'$\alpha=1/2$', r'$\alpha=1/4$', r'$\alpha=1/8$', r'$\alpha=1/16$', r'$\alpha=1/32$']

    plt.switch_backend('pgf')

    fig, ax = plt.subplots(figsize=(4.0, 2.5))
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": "Times New Roman",
    })

    major_ticks_x = [0.1, 1, 10, 100, 1000]
    ax.set_xticks(major_ticks_x)
    minor_ticks_x = np.concatenate([np.linspace(2, 9, 8) * 10**exp for exp in range(-1, 4)])
    ax.xaxis.set_minor_locator(FixedLocator(minor_ticks_x))
    ax.xaxis.set_minor_formatter(NullFormatter())

    major_ticks_y = [1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 1]
    ax.set_yticks(major_ticks_y)
    minor_ticks_y = np.concatenate([np.linspace(2, 9, 8) * 10**exp for exp in range(-8, 0)])
    ax.yaxis.set_minor_locator(FixedLocator(minor_ticks_y))
    ax.yaxis.set_minor_formatter(NullFormatter())

    for j in range(len(Alpha_values)):
        alpha = Alpha_values[j]

        TTime_values = data[f"{alpha}"]["tot_time"]
        Density_values = data[f"{alpha}"]["density"]

        plt.plot(TTime_values, Density_values, marker=marker[j], markersize=markersize[j], markeredgewidth=markeredgewidth[j], label=labels[j], alpha=1.0, linewidth=0.35)

    plt.tick_params(axis='both', which='major', direction='in', labelsize=8, length=5, width=0.2, bottom=True, top=True, left=True, right=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=2.2, width=0.1, bottom=True, top=True, left=True, right=True)

    plt.axis((2.6, 1.2e03, 1.7e-04, 1.6e-01))

    plt.xlabel(r'\bf{total annealing time} \textbf{$\tau_{total}$}', fontsize=10)
    plt.ylabel(r'\bf{excitation density} \textbf{$d$}', fontsize=10)

    legend = plt.legend(loc='lower left', fontsize=7.0, frameon=False, bbox_to_anchor=(0.02, 0.02))

    plt.tight_layout(pad=0.2)
    plt.subplots_adjust(right=0.9)
    plt.savefig(f'plots/inhomo_test_log_plot_{N}.pdf')

    plt.close()


if __name__ == '__main__':

    Alpha_values = [0, 1/2, 1/4, 1/8, 1/16, 1/32] # 0 stands for homogeneous transition
    N = 1000

    filename = f"data/inhomo_test_log_data_{N}.json"

    plot_inhomo_log(N, Alpha_values, filename)