import numpy as np
import matplotlib.pyplot as plt
import json

def plot_inhomo_lin(N, Alpha_values, filename):

    with open(filename, 'r') as f:

            data = json.load(f)

    marker = [r'$\ast$', '+', 'x', r'$\star$', r'$\circ$']
    markersize = [4, 4, 3.5, 4, 3]
    markeredgewidth = [0.25, 0.3, 0.3, 0.2, 0.05]
    labels = [r'$\alpha=1/2$', r'$\alpha=1/4$', r'$\alpha=1/8$', r'$\alpha=1/16$', r'$\alpha=1/32$']

    plt.switch_backend('pgf')

    fig, ax = plt.subplots(figsize=(4.0, 2.5))  # Two-column width in inches, aspect ratio adjusted for better fit

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": "Times New Roman",
    })

    plt.plot([], [], alpha=0.0)  # Dummy plot for consistency with legend

    for j in range(len(Alpha_values)):
        alpha = Alpha_values[j]

        Velocity_values = data[f"{alpha}"]["velocity"]
        Density_values = data[f"{alpha}"]["density"]

        plt.plot(Velocity_values[::2], Density_values[::2], marker=marker[j], markersize=markersize[j], markeredgewidth=markeredgewidth[j], label=labels[j], alpha=1.0, linewidth=0.35)

    major_ticks_y = [0.00, 0.02, 0.04, 0.06, 0.08]
    ax.set_yticks(major_ticks_y)
    y_tick_labels = [r'$0.00$', r'$0.02$', r'$0.04$', r'$0.06$', r'$0.08$']
    ax.set_yticklabels(y_tick_labels)

    major_ticks_x = [0, 2, 4, 6, 8, 10]
    ax.set_xticks(major_ticks_x)
    x_tick_labels = [r'$0$', r'$2$', r'$4$', r'$6$', r'$8$', r'$10$']
    ax.set_xticklabels(x_tick_labels)

    plt.tick_params(axis='both', which='major', direction='in', labelsize=8, length=5, width=0.1, bottom=True, top=True, left=True, right=True)

    x_values = [2.0] * 100
    y_values = np.linspace(0.00, 0.012, 100)
    plt.plot(x_values, y_values, linestyle='--', color='navy', label=r'$c=2$', linewidth=0.7)

    plt.plot([], [], alpha=0.0, label=' ')   # Dummy plot
    plt.plot([], [], alpha=0.0, label=' ')   # Dummy plot

    plt.xlim(0, 10.2)
    plt.ylim(0, 0.08)

    plt.xlabel(r'\bf{spatial ramp velocity} \textbf{$v$}', fontsize=10)
    plt.ylabel(r'\bf{excitation density} \textbf{$d$}', fontsize=10)

    plt.legend(loc='upper left', fontsize=7.7, frameon=False, bbox_to_anchor=(0.02, 0.98), ncol=2)

    plt.tight_layout(pad=0.2)
    plt.subplots_adjust(right=0.9)
    plt.savefig(f'plots/inhomo_test_lin_plot_{N}.pdf')

    plt.close()

if __name__ == '__main__':

    Alpha_values = [1/2, 1/4, 1/8, 1/16, 1/32] # 0 stands for homogeneous transition
    N = 1000

    filename = f"data/inhomo_test_lin_data_{N}.json"

    plot_inhomo_lin(N, Alpha_values, filename)