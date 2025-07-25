import random
import NetfluxODE
import numpy as np
from scipy.integrate import ode, odeint
import InputParameters
import pandas as pd
import os
from openpyxl import load_workbook

def RunNetworkModel(tspan, y0, params, speciesNames, specific_time):
    # Step 1: Run the simulation
    t = []
    dt = tspan[1] / ((specific_time * 150) / 10)
    r = ode(NetfluxODE.ODEfunc).set_integrator('vode', method='adams', order=10, rtol=0, atol=1e-6,
                                               with_jacobian=False)
    r.set_initial_value(y0, tspan[0]).set_f_params(*params)
    results = np.empty([0, len(speciesNames)])

    while r.successful() and r.t <= tspan[1]:
        r.integrate(r.t + dt)
        results = np.append(results, [r.y], axis=0)
        t.append(r.t)

    # Step 2: Extract the results at the specific time
    index = int((specific_time - tspan[0]) / (tspan[1] / ((specific_time * 150) / 10))) - 1
    results_at_specific_time = results[index, :]

    # Step 3: Create and return a dictionary of results
    results_dict = {speciesNames[i]: results_at_specific_time[i] for i in range(len(speciesNames))}
    return results_dict


def RunABMModel(results_dict, TIME_STEPS, w_IL1B, w_IL6, w_TNFA, collagen):
    def ode_system_latentTGFB(latentTGFB_ABM, t, kgen, ksec, latentTGFBnet, kdeg, kact):
        return kgen + ksec * latentTGFBnet - kdeg * latentTGFB_ABM - kact * latentTGFB_ABM

    # Set constants
    kgen, ksec, latentTGFBnet, kdeg, kact = 530000, 23700, results_dict['latentTGFB'], 0.0096, 0.045

    # Set initial condition and time points
    initial_latentTGFBABM, t = results_dict["latentTGFB"], np.linspace(0, TIME_STEPS)

    # Solve the ODE
    latentTGFB_ABM = odeint(ode_system_latentTGFB, initial_latentTGFBABM, t, args=(kgen, ksec, latentTGFBnet, kdeg, kact))[-1][0]
    ########################################################################################################################
    # Equation Solution TGFB_ABM (6)
    TGFB_ABM = 0.045 * latentTGFB_ABM

    ########################################################################################################################
    # Differential Equation Solution IL-1B (7)
    def ode_system_IL1B(IL1B_ABM, t, kgen, kdeg):
        return kgen - kdeg * IL1B_ABM

    # Set constants
    kgen, kdeg = 4847, 0.277

    # Set initial condition and time points
    initial_IL1B, t = w_IL1B, np.linspace(0, TIME_STEPS)

    # Solve the ODE
    IL1b_ABM = odeint(ode_system_IL1B, initial_IL1B, t, args=(kgen, kdeg))[-1][0]

    ########################################################################################################################
    # Differential Equation Solution IL-6_ABM (8)
    def ode_system_IL6(IL6_ABM, t, kgen, ksec, IL6net, kdeg):
        return kgen + ksec * IL6net - kdeg * IL6_ABM

    # Set constants
    kgen, ksec, kdeg, IL6net = 256000, 79360, 0.277, results_dict['IL6']

    # Set initial condition and time points
    initial_IL6, t = w_IL6, np.linspace(0, TIME_STEPS)

    # Solve the ODE
    IL6_ABM = odeint(ode_system_IL6, initial_IL6, t, args=(kgen, ksec, IL6net, kdeg))[-1][0]

    ########################################################################################################################
    # Differential Equation Solution TNFA_ABM (9)
    def ode_system_TNFA(TNFA_ABM, t, kgen, kdeg):
        return kgen - kdeg * TNFA_ABM

    # Set constants
    kgen, kdeg = 895.4, 1.386

    # Set initial condition and time points
    initial_TNFA, t = w_TNFA, np.linspace(0, TIME_STEPS)

    # Solve the ODE
    TNFA_ABM = odeint(ode_system_TNFA, initial_TNFA, t, args=(kgen, kdeg))[-1][0]

    ########################################################################################################################
    # Differential Equation Solution Collagen (10)
    def ode_system_Collagen(Collagen_ABM, t, kdep, kdeg, ColIRNAnet, ColIIIRNAnet):
        return kdep * (ColIRNAnet + ColIIIRNAnet) - kdeg * Collagen_ABM

    # Set constants
    kdep, kdeg, ColIRNAnet, ColIIIRNAnet = 0.0056, 0.0035, results_dict['CImRNA'], results_dict['CIIImRNA']

    # Set initial condition and time points
    initial_Collagen, t = collagen, np.linspace(0, TIME_STEPS)

    # Solve the ODE
    Collagen_ABM = odeint(ode_system_Collagen, initial_Collagen, t, args=(kdep, kdeg, ColIRNAnet, ColIIIRNAnet))[-1][0]

    w_IL6 = IL6_ABM / (IL6_ABM + 462000)  # kd_IL6 = 462000
    w_IL1B = IL1b_ABM / (IL1b_ABM + 8750)  # kd_IL1B = 8750
    w_TNFA = TNFA_ABM / (TNFA_ABM + 323)  # kd_TNFA = 32
    w_TGFB = TGFB_ABM / (TGFB_ABM + 700)  # kd_TGFB = 700

    return [w_IL6, w_IL1B, w_TNFA, w_TGFB, Collagen_ABM]


def RunSimulation(w, TIME_STEPS, SPECIFIC_TIME, num):
    # Set Placeholders for ABM Equation Variables
    w_TGFB = w_IL6 = w_IL1B = w_TNFA = 0
    collagen = 0
    results_dict = {}

    for x in range(TIME_STEPS):
        if x != 0:
            w[1] = w_TGFB  # TGFB
            w[3] = w_IL6  # IL6
            w[4] = w_IL1B  # IL1
            w[5] = w_TNFA  # TNFa
            y0 = list(results_dict.values())

        results_dict = RunNetworkModel([0, SPECIFIC_TIME], InputParameters.y0,
                                       [InputParameters.tau, InputParameters.ymax,
                                        w, InputParameters.n, InputParameters.EC50],
                                       InputParameters.speciesNames, SPECIFIC_TIME)

        resulting_parameters = RunABMModel(results_dict, TIME_STEPS, w_IL1B, w_IL6, w_TNFA, collagen)

        w_TGFB = resulting_parameters[3]
        w_IL6 = resulting_parameters[0]
        w_IL1B = resulting_parameters[1]
        w_TNFA = resulting_parameters[2]
        collagen = resulting_parameters[4]

    if num > 0:
        print("Simulation " + str(num) + " Collagen: " + str(collagen))
    else:
        print(f'Actual Collagen with Predicted Values: {collagen}')
    return collagen

def generate_simulation_data(num_simulations, time_steps, specific_time, w_base, file_name="simulation_data.xlsx"):
    # Base input structure for simulations
    w_base = [0.25] * 11 + [1] * (len(w_base) - 11)  # Adjust the length based on your requirements

    # Check if file exists to set up headers if needed
    file_exists = os.path.exists(file_name)

    for sim in range(num_simulations):
        # Randomize the first 11 entries in the input parameters (w)
        w = w_base.copy()
        for i in range(11):
            w[i] = random.uniform(0, 1)

        # Run simulation and store results
        collagen = RunSimulation(w, time_steps, specific_time, sim + 1)

        data_row = {
            'AngII': w[0], 'TGFB': w[1], 'mechanical': w[2], 'IL6': w[3], 'IL1': w[4],
            'TNFa': w[5], 'NE': w[6], 'PDGF': w[7], 'ET1': w[8], 'NP': w[9],
            'Forskolin': w[10], 'Collagen': collagen
        }
        data_df = pd.DataFrame([data_row])

        if file_exists:
            with pd.ExcelWriter(file_name, engine="openpyxl", mode="a", if_sheet_exists="overlay") as writer:
                data_df.to_excel(writer, index=False, header=False, startrow=writer.sheets['Sheet1'].max_row)
        else:
            # Write with headers if file does not exist
            data_df.to_excel(file_name, index=False)
            file_exists = True  # Set flag so future writes will append without headers

        print(f"Simulation {sim + 1} completed and saved to {file_name}")
        print("===================================================================")

def main():
    w_base = [2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01,
         2.500000e-01, 2.500000e-01, 2.500000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]

    TIME_STEPS = 10
    SPECIFIC_TIME = 10

    generate_simulation_data(num_simulations=10000, time_steps=TIME_STEPS,
                             specific_time=SPECIFIC_TIME, w_base=w_base)

if __name__ == "__main__":
    main()
