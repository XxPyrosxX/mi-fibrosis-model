from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.model_selection import train_test_split
import sys
sys.path.append('../model')
from RunSimulations import RunSimulation
# Uses regular random forrest to generate feature importance
# It isn't as efficient in sampling the best inputs

def run_simulation_with_variation(param_index, start, end, step, w_base):

    end += 0.05

    # Initialize base parameters
    TGFB = 0.50
    mechanical = 0.70
    IL6 = 0.60
    IL1 = 0.25
    TNFa = 0.35
    NE = 0.85
    PDGF = 0.92
    ET1 = 0.07
    NP = 0.15
    Forskolin = 0.15

    param_indices = {
    0: 'AngII',
    1: 'TGFB',
    2: 'mechanical',
    3: 'IL6',
    4: 'IL1',
    5: 'TNFa',
    6: 'NE',
    7: 'PDGF',
    8: 'ET1',
    9: 'NP',
    10: 'Forskolin'
    }

    # List of base parameters
    base_params = [0, TGFB, mechanical, IL6, IL1, TNFa, NE, PDGF, ET1, NP, Forskolin]

    # Generate the range of values for the specified parameter
    param_values = [round(val, 2) for val in frange(start, end, step)]
    param_name = param_indices[param_index]
    # Run simulations with varying parameter
    for value in param_values:
        w = w_base.copy()
        w[0:11] = base_params.copy()  # Ensure all parameters are set to base values
        w[param_index] = value  # Vary the specified parameter

        collagen = RunSimulation(w, 10, 10, 1)
        print(f"Parameter {param_name}: {value:.2f}, Collagen: {collagen}")
    print("---------------------------------------")

# Helper function to generate float ranges
def frange(start, end, step):
    while start <= end:
        yield start
        start += step

def main():

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

    file_path = '../data/trial20.xlsx'
    data = pd.read_excel(file_path)
    X = data.drop(columns=["Collagen"])
    y = data[["Collagen"]]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    rf_model = RandomForestRegressor(n_estimators=100, random_state=42)
    rf_model.fit(X_train, y_train)

    feature_importances = rf_model.feature_importances_

    feature_importance_df = pd.DataFrame({'Feature': X.columns, 'Importance': feature_importances})
    feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=False)
    print(feature_importance_df)

    # Step 7: Plot feature importance - Upgraded Visualization
    plt.figure(figsize=(18, 12))
    sns.set(style="whitegrid", font_scale=1.4, rc={
        "axes.edgecolor": "black",
        "axes.linewidth": 1.25,
        "grid.linewidth": 0.5,
        "axes.labelweight": "bold",
        "axes.titlesize": 20,
        "axes.titleweight": "bold"
    })

    ax = sns.barplot(
        x='Importance',
        y='Feature',
        data=feature_importance_df,
        palette='viridis'
    )

    # Add importance values on each bar
    for i, (importance, feature) in enumerate(
            zip(feature_importance_df["Importance"], feature_importance_df["Feature"])):
        ax.text(importance + 0.002, i, f'{importance:.3f}', va='center', ha='left', fontsize=10, color='black', weight='semibold')

    plt.title('Feature Importance for Collagen Production', fontsize=16, weight='bold', pad=15)
    plt.xlabel('Importance Score', fontsize=13, labelpad=10, weight="bold")
    plt.ylabel('Parameter', fontsize=13, labelpad=10, weight="bold")
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.tight_layout()
    sns.despine()
    plt.show()


    print("---------------------------------------")
    index = int(input("Enter Index of Param you want to Change: "))
    run_simulation_with_variation(param_index=index, start=0, end=1.0, step=0.05, w_base=w_base)
    print("---------------------------------------")


if __name__ == "__main__":
    main()