import matplotlib.pyplot as plt
import seaborn as sns
from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import pandas as pd
import sys
sys.path.append('../model')
from RunSimulations import RunSimulation

# Setup Seaborn style for aesthetics
sns.set(style="whitegrid")

# TOO MANY SIMULATIONS NEEDED TO RUN THIS METHOD ACCURATELY
# Global Sensitivity Analysis (GSA): Evaluates how variations in all parameters affect output across the entire input space.
# Sobol Analysis:
#   function generates N x (2D + 2) samples where D is # of parameters and N is base sample size
#   Ran with N = 2, and got fairly accurate results (48 runs of simulation)
#   Similar to Random Forrest but has a better sampling of the input instead of relying on a dataset
# Morris Analysis: Takes way too many inputs to get accurate results

def run_sensitivity_analysis():

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

    problem = {
        'num_vars': 11,
        'names': ['AngII', 'TGFB', 'mechanical', 'IL6', 'IL1', 'TNFa', 'NE', 'PDGF', 'ET1', 'NP', 'Forskolin'],
        'bounds': [[0, 1]] * 11  # Assuming all parameters vary between 0 and 1
    }

    # Generate samples for sensitivity analysis
    param_values = saltelli.sample(problem, 2)

    # Run the simulation for each sample
    collagen_outputs = np.array(
        [RunSimulation(list(params) + [1] * (len(w_base) - 11), 10, 10, 1) for params in param_values])

    # Sobol analysis
    sobol_results = sobol.analyze(problem, collagen_outputs)
    plot_sobol_results(sobol_results, problem['names'])

    # # Morris method
    # param_values_morris = morris.sample(problem, 1024, num_levels=4)
    # collagen_outputs_morris = np.array(
    #     [RunSimulation(list(params) + [1] * (len(w_base) - 11), 10, 10, 1) for params in param_values_morris])
    # morris_results = morris_analyze.analyze(problem, param_values_morris, collagen_outputs_morris)
    # print("Morris Sensitivity Analysis Results:")
    # print("Mean:", morris_results['mu'])
    # print("Mean Absolute Deviation:", morris_results['mu_star'])
    # print("Standard Deviation:", morris_results['sigma'])


def plot_sobol_results(sobol_results, parameter_names):
    """
    Visualize Sobol Sensitivity Analysis results using bar plots for both first and total order sensitivity indices.
    Show total order values above each bar.
    """
    # Prepare DataFrame for First Order Sensitivity (S1)
    first_order_df = pd.DataFrame({
        'Parameter': parameter_names,
        'First Order Sensitivity (S1)': sobol_results['S1'],
    })

    # Prepare DataFrame for Total Order Sensitivity (ST)
    total_order_df = pd.DataFrame({
        'Parameter': parameter_names,
        'Total Order Sensitivity (ST)': sobol_results['ST'],
    })

    # Plot First Order Sensitivity Indices (S1)
    plt.figure(figsize=(10, 6))
    ax1 = sns.barplot(data=first_order_df, x='Parameter', y='First Order Sensitivity (S1)', palette="Blues_d")
    plt.title('First Order Sensitivity Indices (S1)', fontsize=16)
    plt.xlabel('Parameter', fontsize=12)
    plt.ylabel('S1 Value', fontsize=12)
    plt.xticks(rotation=45, ha='right')

    # Annotate each bar with the first order value
    for p in ax1.patches:
        ax1.annotate(f'{p.get_height():.4f}',
                     (p.get_x() + p.get_width() / 2., p.get_height()),
                     ha='center', va='center',
                     fontsize=10, color='black',
                     xytext=(0, 8), textcoords='offset points')

    plt.tight_layout()
    plt.show()

    # Plot Total Order Sensitivity Indices (ST) with values displayed above each bar
    plt.figure(figsize=(10, 6))
    ax2 = sns.barplot(data=total_order_df, x='Parameter', y='Total Order Sensitivity (ST)', palette="Oranges_d")
    plt.title('Total Order Sensitivity Indices (ST)', fontsize=16)
    plt.xlabel('Parameter', fontsize=12)
    plt.ylabel('ST Value', fontsize=12)
    plt.xticks(rotation=45, ha='right')

    # Annotate each bar with the total order value
    for p in ax2.patches:
        ax2.annotate(f'{p.get_height():.4f}',
                     (p.get_x() + p.get_width() / 2., p.get_height()),
                     ha='center', va='center',
                     fontsize=10, color='black',
                     xytext=(0, 8), textcoords='offset points')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_sensitivity_analysis()