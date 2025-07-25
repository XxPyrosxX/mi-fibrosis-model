import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.append('../model')
from RunSimulations import RunSimulation

def monte_carlo_simulation(num_samples=1000):
    """
    Perform Monte Carlo simulation by sampling parameter values from probability distributions
    and evaluating their effect on collagen output.
    """

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

    # Define mean and standard deviation for each parameter (assume normal distribution)
    param_means = np.array([0.50, 0.70, 0.60, 0.25, 0.35, 0.85, 0.92, 0.07, 0.15, 0.15, 0.20])
    param_stds = np.array([0.1] * len(param_means))  # Assuming 10% standard deviation

    # Generate samples from normal distributions
    sampled_params = np.random.normal(loc=param_means, scale=param_stds, size=(num_samples, len(param_means)))

    collagen_outputs = []

    for i in range(num_samples):
        w_base[:11] = sampled_params[i]
        collagen = RunSimulation(w_base, 10, 10, 1)
        collagen_outputs.append(collagen)

    # Convert to DataFrame for analysis
    df_results = pd.DataFrame(sampled_params,
                              columns=['TGFB', 'Mechanical', 'IL6', 'IL1', 'TNFa', 'NE', 'PDGF', 'ET1', 'NP',
                                       'Forskolin', 'AngII'])
    df_results['Collagen'] = collagen_outputs

    return df_results


def analyze_uncertainty(df_results):
    """Perform statistical analysis and visualization of uncertainty in collagen predictions."""
    print("Summary statistics:")
    print(df_results.describe())

    # Pairplot to visualize relationships (optional, can comment this out for cleaner presentations)
    # sns.pairplot(df_results, diag_kind='kde')
    # plt.show()

    # Professional Correlation Heatmap
    plt.figure(figsize=(18, 12))
    corr = df_results.corr()

    # Mask to show only upper triangle
    mask = np.triu(np.ones_like(corr, dtype=bool))

    # Custom diverging colormap
    cmap = sns.diverging_palette(220, 20, as_cmap=True)

    sns.heatmap(
        corr,
        mask=mask,
        cmap=cmap,
        vmax=1.0,
        vmin=-1.0,
        center=0,
        annot=True,
        fmt=".2f",
        square=True,
        linewidths=.5,
        cbar_kws={"shrink": 0.8, "label": "Correlation Coefficient"},
        annot_kws={"size": 10}
    )

    plt.title('Correlation Between Parameters and Collagen Output', fontsize=18, weight='bold', pad=20)
    plt.xticks(fontsize=10, rotation=45, ha='right')
    plt.yticks(fontsize=10)
    plt.tight_layout()
    plt.show()


def main():
    num_samples = 25  # Define number of Monte Carlo samples
    df_results = monte_carlo_simulation(num_samples=num_samples)
    analyze_uncertainty(df_results)


if __name__ == "__main__":
    main()
