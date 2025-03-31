from ising_model.ising_model import Ising
import json
import multiprocessing

g_f = 1.0
g_i = 2.0

dt = 0.01

L = 1000

Delta_TTime = 4/3
TTime_values = [Delta_TTime ** i for i in range(-10, 25)]
Alpha_values = [0, 1/2, 1/4, 1/8, 1/16, 1/32] # alpha=0 represents homogeneous transition

system = Ising(L, g_i, g_f)

def compute_density_for_alpha(alpha):

    if alpha==0:
        # homogeneous phase transition
        Density_values = []

        for TTime in TTime_values:
            print(f'total_time = {TTime}')
            Density_values.append(system.homo_evolution(TTime, dt))

    else:
        # inhomogeneous phase transition
        Density_values = []

        for TTime in TTime_values:
            print(f'total_time = {TTime}')
            Density_values.append(system.inhomo_evolution(TTime, dt, alpha))

    return alpha, {"tot_time": TTime_values, "density": Density_values}

if __name__ == '__main__':

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(compute_density_for_alpha, Alpha_values)

    # save data to a file

    data = dict(results)
    filename = f"data/inhomo_test_log_data_{L}.json"

    with open(filename, 'w') as f:
        json.dump(data, f)
