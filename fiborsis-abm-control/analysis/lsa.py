from SALib.analyze import sobol
from SALib.sample import saltelli
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('../model')
from RunSimulations import RunSimulation

# Most sensitive impacts to the collagen
# The parameters with the most effect even with a very small change in the parameter

w_base = [2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01,
          2.500000e-01,
          2.500000e-01, 2.500000e-01, 2.500000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1,
          1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1,
          1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
          1, 1,
          1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]

def collagen_simulation(w):
    """ Placeholder function to simulate collagen production. """
    return RunSimulation(w, 10, 10, 1)


def run_sobol_analysis():

    problem = {
        'num_vars': 11,
        'names': ['AngII', 'TGFB', 'mechanical', 'IL6', 'IL1', 'TNFa', 'NE', 'PDGF', 'ET1', 'NP', 'Forskolin'],
        'bounds': [[0, 1]] * 11
    }

    # Generate samples
    param_values = saltelli.sample(problem, 2)

    # Run model
    Y = np.array([RunSimulation(list(params) + [1] * (len(w_base) - 11), 10, 10, 1) for params in param_values])

    # Perform Sobol analysis
    sobol_indices = sobol.analyze(problem, Y)

    # Print and plot results
    print("Sobol First Order Indices:", sobol_indices['S1'])
    print("Sobol Total Order Indices:", sobol_indices['ST'])

    sns.barplot(x=problem['names'], y=sobol_indices['ST'])
    plt.xticks(rotation=45)
    plt.title("Sobol Sensitivity Analysis (Total Order)")
    plt.show()


def run_local_sensitivity_analysis():
    problem = {
        'num_vars': 11,
        'names': ['AngII', 'TGFB', 'mechanical', 'IL6', 'IL1', 'TNFa', 'NE', 'PDGF', 'ET1', 'NP', 'Forskolin'],
        'bounds': [[0, 1]] * 11
    }

    w_new = w_base.copy()
    w_new[:11] = [0.5] * 11
    base_point = np.array(w_new)  # Central point in parameter space
    delta = 0.01  # Small perturbation

    sensitivities = []

    base_output = collagen_simulation(base_point)

    for i in range(problem['num_vars']):
        perturbed_point = base_point.copy()
        perturbed_point[i] += delta
        perturbed_output = collagen_simulation(perturbed_point)

        local_sensitivity = (perturbed_output - base_output) / delta
        sensitivities.append(local_sensitivity)

    sns.barplot(x=problem['names'], y=sensitivities)
    plt.xticks(rotation=45)
    plt.title("Local Sensitivity Analysis")
    plt.show()


if __name__ == "__main__":
    #run_sobol_analysis()
    run_local_sensitivity_analysis()
