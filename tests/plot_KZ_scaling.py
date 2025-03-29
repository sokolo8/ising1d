import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter
from sklearn.linear_model import LinearRegression
from scipy import stats
import json


def plot_KZ_scaling(N_values, filename):

    with open(filename, 'r') as f:

        data = json.load(f)

    marker = [r'$\otimes$', r'$\ast$', '+', 'x', r'$\star$', r'$\circ$']
    markersize = [3, 4, 4, 3.5, 4, 3]
    markeredgewidth = [0.25, 0.25, 0.3, 0.3, 0.2, 0.05]

    plt.switch_backend('pgf')

    # Main plot
    fig, ax = plt.subplots(figsize=(4.0, 2.5))
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": "Times New Roman",
    })

    major_ticks_x = [0.01, 0.1, 1, 10, 100, 1e03, 1e04]
    ax.set_xticks(major_ticks_x)
    minor_ticks_x = np.concatenate([np.linspace(2, 9, 4) * 10**exp for exp in range(-2, 5)])
    ax.xaxis.set_minor_locator(FixedLocator(minor_ticks_x))
    ax.xaxis.set_minor_formatter(NullFormatter())

    major_ticks_y = [1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01]
    ax.set_yticks(major_ticks_y)
    minor_ticks_y = np.concatenate([np.linspace(2, 8, 4) * 10**exp for exp in range(-8, 0)])
    ax.yaxis.set_minor_locator(FixedLocator(minor_ticks_y))
    ax.yaxis.set_minor_formatter(NullFormatter())

    plt.tick_params(axis='both', which='major', direction='in', labelsize=8, length=5, width=0.2, bottom=True, top=True, left=True, right=True)
    plt.tick_params(axis='both', which='minor', direction='in', length=2.2, width=0.1, bottom=True, top=True, left=True, right=True)

    # Loop for main plot
    for j in range(len(N_values)):
        N = N_values[j]

        TTime_values = data[f"{N}"]["tot_time"]
        Density_values = data[f"{N}"]["density"]

        plt.plot(TTime_values, Density_values, marker=marker[j], markersize=markersize[j], markeredgewidth=markeredgewidth[j], label=rf'$N = {N}$', alpha=1.0, linewidth=0.35)

    # Linear regression for max N

    T_min = 11
    T_max = 26

    TTime_values_log = np.log10(TTime_values)
    Density_values_log = np.log10(Density_values)

    Density_values_log = Density_values_log[T_min:T_max]
    TTime_values_log = TTime_values_log[T_min:T_max]

    TTime_values_log = np.array(TTime_values_log).reshape((-1, 1))

    regression = LinearRegression().fit(np.real(TTime_values_log), np.real(Density_values_log))

    TTime_values_log = TTime_values_log.flatten()

    a = regression.coef_[0]
    b = regression.intercept_

    # Calculating slope error

    residuals = np.real(Density_values_log) - b - a * np.real(TTime_values_log)

    n = len(np.real(TTime_values_log))
    p = 2
    df = n - p

    residual_variance = np.sum(residuals ** 2) / df

    SE_beta = np.sqrt(residual_variance / np.sum((np.real(TTime_values_log) - np.mean(np.real(TTime_values_log))) ** 2))

    confidence_level = 0.95
    t_value = np.abs(stats.t.ppf((1 - confidence_level) / 2, df))

    confidence_interval = t_value * SE_beta

    t_min = min(TTime_values_log)
    t_max = max(TTime_values_log)
    x_values = np.linspace(t_min, t_max, 100)
    y_values = a * x_values + b
    x_values = 10 ** x_values
    y_values = 10 ** y_values
    plt.plot(x_values, y_values, linestyle='--', color='navy', label=rf'lin fit ${max(N_values)}$, $a={a:.2f}({confidence_interval * 100 :.0f})$', linewidth=1.3)


    plt.axis((0.035, 1.3e03, 1e-06, 2*1e-01))
    plt.xlabel(r'\bf{total annealing time} $\tau_{total}$', fontsize=10)
    plt.ylabel(r'\bf{excitation density} $d$', fontsize=10)
    plt.legend(loc='lower left', fontsize=7.0, frameon=False, bbox_to_anchor=(0.02, 0.02))


    plt.tight_layout(pad=0.2)
    plt.subplots_adjust(right=0.9)
    plt.savefig('plots/KZ_scaling_plot.pdf')

    plt.close()


if __name__ == '__main__':

    N_values = [100, 200, 300, 400, 500, 1000]

    filename = "data/KZ_scaling_data.json"

    plot_KZ_scaling(N_values, filename)