from ising_model import Ising
import json
import multiprocessing

g_f = 1.0
g_i = 2.0

dt = 0.01

N_values = [100, 200, 300, 400, 500, 1000]

Delta_TTime = 4/3
TTime_values = [Delta_TTime ** i for i in range(-10, 25)]

# homogeneous phase transition

def compute_density_for_N(N):
    system = Ising(N, g_i, g_f)
    Density_values = []

    for TTime in TTime_values:
        print(f'total_time = {TTime}')
        Density_values.append(system.homo_evolution(TTime, dt))

    return N, {"tot_time": TTime_values, "density": Density_values}

if __name__ == "__main__":

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(compute_density_for_N, N_values)

    # save data to a file

    data = dict(results)
    filename = "data/KZ_scaling_data.json"

    with open(filename, 'w') as f:
        json.dump(data, f)
