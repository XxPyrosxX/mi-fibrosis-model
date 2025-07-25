# 🧬 Control Methods for Dynamic Cellular Simulations of Fibrosis

**Rohun Gargya**  
College of Medicine, University of Florida  
**Advisor:** Dr. Luis Sordo Vieira, Department of Medicine

---

## 📘 Overview

This repository contains the codebase, simulations, and analysis pipeline for developing control and sensitivity methods for dynamic cellular simulations of fibrosis. We focus on enhancing the efficiency and interpretability of Agent-Based Models (ABMs) coupled with Ordinary Differential Equations (ODEs) by integrating machine learning and global sensitivity techniques.

Our approach enables accurate parameter importance ranking, uncertainty quantification, and prediction of fibrosis outcomes (e.g., collagen production) with significantly reduced simulation costs.

---

## 🧪 Abstract

Agent-based models are powerful tools for modeling dynamic biological systems, but their high computational cost and complexity can limit usability in large-scale analysis. This project presents a hybrid framework integrating:

- 🌲 **Random Forest Regression**
- 📊 **Sobol Sensitivity Analysis**
- 🧠 **SHAP (Shapley Additive Explanations)**
- 🎲 **Monte Carlo Simulation**

...to approximate and analyze fibrosis simulations while minimizing the number of required runs. We apply our methods to a multi-scale ABM-ODE model of fibrosis, identifying key parameter drivers of collagen production, a central process in fibrotic progression.

---

## 🔬 Model Description

- **ABM Component:** Models extracellular cellular behavior and cytokine signaling.
- **ODE Component:** Intracellular signaling network with 91 nodes and 142 reactions.
- **Inputs:** IL-6, IL-1β, TNFα, TGFβ (as pg/mL).
- **Outputs:** Collagen mRNA, protein markers, normalized activity weights (0–1).
- **Simulation Period:** Simulates a 6-week time course.

Original MATLAB models were ported to **Python** to take advantage of open-source tools and reduce runtime complexity.

---

## 🧰 Analysis Methods

| Method | Description | Purpose |
|--------|-------------|---------|
| **Random Forest Regressor** | ML model to predict collagen levels | Fast prediction with limited training data |
| **Sobol Sensitivity Analysis** | Variance-based method for parameter ranking | Quantifies first- and total-order effects |
| **SHAP Values** | Model interpretability via Shapley values | Local + global feature impact analysis |
| **Monte Carlo Simulation** | Random sampling of parameter space | Uncertainty estimation and distribution analysis |

These methods are complementary and together provide a full diagnostic of model behavior.

---

## 📁 Project Structure

fibrosis-abm-control/
├── data/                    
│   ├── inputs/              # Input cytokine parameters for simulations
│   └── outputs/             # Simulation results and intermediate data
│
├── model/                   
│   ├── run_simulation.py    # Executes the ABM-ODE simulation pipeline
│   ├── abm_module.py        # Agent-based model logic
│   └── ode_module.py        # ODE signaling network implementation
│
├── analysis/                
│   ├── random_forest.py     # Trains and evaluates the Random Forest Regressor
│   ├── sobol_analysis.py    # Computes Sobol sensitivity indices using SALib
│   ├── shap_explainer.py    # Visualizes feature importance using SHAP
│   └── monte_carlo.py       # Runs Monte Carlo sampling for uncertainty analysis
│
├── notebooks/               
│   ├── exploratory_analysis.ipynb   # Data visualization and inspection
│   ├── shap_visuals.ipynb          # Interactive SHAP plots and force plots
│   └── monte_carlo_results.ipynb   # Output summaries from MC simulations
│
├── results/                 
│   ├── figures/             # Plots and visualizations of analysis results
│   └── summaries/           # CSV/JSON reports, ranked feature lists, etc.
│
├── requirements.txt         # All required Python packages and dependencies
└── README.md                # Overview, usage, and documentation (this file)
