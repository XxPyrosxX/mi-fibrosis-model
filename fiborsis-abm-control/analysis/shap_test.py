import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
import pandas as pd
import shap

# Basic feature to feature comparison using SHAP (Shapley Additive Explanations)
# a method used to explain the output of machine learning models by assigning importance values to each feature.
# It is based on Shapley values, a concept from cooperative game theory that fairly distributes contributions among
# players (features, in this case)

# Load your dataset
file_path = "../data/trial50.xlsx"
data = pd.read_excel(file_path)

X = data.drop(columns=["Collagen"])
y = data["Collagen"]

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a model
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# SHAP Explainer
# The graph displays both positive and negative SHAP interaction values, depending on the combination parameters
# A positive value means the interaction between these features positively affects the model output, while a negative
# value indicates a negative influence.
explainer = shap.TreeExplainer(model)
shap_values = explainer(X_test)

# Compute Interaction Values
interaction_values = explainer.shap_interaction_values(X_test)

# Select two features to analyze (change as needed)
feature_1 = "TGFB"
feature_2 = "IL6"

# Get feature indices
idx1 = list(X_test.columns).index(feature_1)
idx2 = list(X_test.columns).index(feature_2)

# Plot interaction effects
shap.dependence_plot((idx1, idx2), interaction_values, X_test)
plt.show()
