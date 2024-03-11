import numpy as np
import math as m
import matplotlib.pyplot as plt
from tabulate import tabulate
import pandas as pd
import pprint
import warnings

from Utilities.ChemicalProperties import ChemicalProperties
from Utilities.Model import Model
from Utilities.utilities import ToMetric, ToEnglish, plot_graph

# Use for Rhode Setup
initialInputs = {
    # Purpose:  Dictionary of Initial Inputs, organized by section
    "T_tank": ToMetric(280, 'K'),
    "m_T": ToMetric(2.1, 'kg'),
    "m_N2O": ToMetric(0.4, 'lbm'),
    #### Oxidizer Tank
    "V_tank": ToMetric(0.259, 'L'),
    "P_tank": ToMetric(1159, 'psi'),
    #### Injector
    "A_inj": ToMetric(2.51390206896e-6, 'm^2'),
    # "A_inj": ToMetric(3.4e-6, 'm^2'),
    "C_d": ToMetric(0.25, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToMetric(2089.83281, 'kg/m^3'),
    "M_chmb": ToMetric(30.169, 'g/mol'),
    "gamma": ToMetric(1.2516, 'unitless'),
    "T_chmb": ToMetric(4000, 'K'),
    ## Fuel Regression Properties
    "n": ToMetric(1.681, 'unitless'),
    "a": ToMetric(9.33E-8, 'unitless'),
    #### Fuel Grain
    "L_fuel": ToMetric(5.5, 'in'),
    "OD_fuel": ToMetric(2, 'in'),
    "Cylindrical": True,
    "ID_fuel": ToMetric(0.5, 'in'),
    #### Nozzle
    "d_t": ToMetric(0.3, 'in'),
    "alpha": np.deg2rad(15), # Nozzle Diverging Half-Cone Angle
    "d_e": ToMetric(0.77, 'in'),
    ### Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

# Aluminum model
initialInputs1 = {
    ## Purpose:  Dictionary of Initial Inputs, organized by section
    "T_tank": ToMetric(284.817, 'K'),
    "m_T": ToMetric(2.1, 'kg'),
    "m_N2O": ToMetric(1.5, 'kg'),
    #### Oxidizer Tank
    "V_tank": ToMetric(3, 'L'),
    "P_tank": ToMetric(3500, 'psi'),
    #### Injector
    "A_inj": ToMetric(8.70966e-6, 'm^2'),
    "C_d": ToMetric(0.4, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToMetric(2089.83281, 'kg/m^3'),
    "M_chmb": ToMetric(30.169, 'g/mol'),
    "gamma": ToMetric(1.2516, 'unitless'),
    "T_chmb": ToMetric(4000, 'K'),
    # "M_chmb": ToMetric(31.071, 'g/mol'),
    # "gamma": ToMetric(1.2025, 'unitless'),
    # "T_chmb": ToMetric(3665.0, 'K'),
    #### Fuel Regression Properties
    "n": ToMetric(1.681, 'unitless'),
    "a": ToMetric(9.33E-8, 'unitless'),
    #### Fuel Grains
    "L_fuel": ToMetric(10.6, 'in'),
    "OD_fuel": ToMetric(3.375, 'in'),
    "Cylindrical": False,
    ## "ID_fuel": ToMetric(2.75, 'in'),
    #### Nozzle
    "d_t": ToMetric(0.889, 'in'),
    "alpha": np.deg2rad(15), # Nozzle Diverging Half-Cone Angle
    "d_e": ToMetric(2.3, 'in'),
    #### Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

# Aluminum model high end prediction
initialInputs = {
    ## Purpose:  Dictionary of Initial Inputs, organized by section
    "T_tank": ToMetric(280, 'K'),
    "m_T": ToMetric(2.1, 'kg'),
    "m_N2O": ToMetric(2.2, 'kg'),
    #### Oxidizer Tank
    "V_tank": ToMetric(3, 'L'),
    "P_tank": ToMetric(3000, 'psi'),
    #### Injector
    "A_inj": ToMetric(1.8e-5, 'm^2'),
    "C_d": ToMetric(0.4, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToMetric(2089.83281, 'kg/m^3'),
    "M_chmb": ToMetric(30.169, 'g/mol'),
    "gamma": ToMetric(1.2516, 'unitless'),
    "T_chmb": ToMetric(4000, 'K'),
    # "M_chmb": ToMetric(31.071, 'g/mol'),
    # "gamma": ToMetric(1.2025, 'unitless'),
    # "T_chmb": ToMetric(3665.0, 'K'),
    #### Fuel Regression Properties
    "n": ToMetric(1.681, 'unitless'),
    "a": ToMetric(9.33E-8, 'unitless'),
    #### Fuel Grains
    "L_fuel": ToMetric(11.5, 'in'),
    "OD_fuel": ToMetric(3.375, 'in'),
    "Cylindrical": False,
    ## "ID_fuel": ToMetric(2.75, 'in'),
    #### Nozzle
    "d_t": ToMetric(0.889, 'in'),
    "alpha": np.deg2rad(15), # Nozzle Diverging Half-Cone Angle
    "d_e": ToMetric(2.3, 'in'),
    #### Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

# Alternate Fuel Testing
initialInputs2 = {
    # Purpose:  Dictionary of Initial Inputs, organized by section
    "T_tank": ToMetric(280, 'K'),
    "m_T": ToMetric(2.1, 'kg'),
    "m_N2O": ToMetric(5, 'lbm'),
    #### Oxidizer Tank
    "V_tank": ToMetric(3, 'L'),
    "P_tank": ToMetric(3000, 'psi'),
    #### Injector
    "A_inj": ToMetric(1.8e-5, 'm^2'),
    "C_d": ToMetric(0.4, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToMetric(920, 'kg/m^3'),
    "M_chmb": ToEnglish(23.719/1000, 'unitless'),
    "gamma": ToEnglish(1.2641, 'unitless'),
    "T_chmb": ToEnglish(3075, 'unitless'),
    ## Fuel Regression Properties
    # "n": ToMetric(1.6386, 'unitless'),
    # "a": ToMetric(9.33E-8, 'unitless'),
    # "n": ToEnglish(0.3667, 'unitless'),
    # "a": ToEnglish(3.413500729493111e-05, 'unitless'),
    # "n": ToEnglish(0.57, 'unitless'),
    # "a": ToEnglish(9.2e-05, 'unitless'),
    "n": ToEnglish(0.67, 'unitless'),
    "a": ToEnglish(0.104/1000, 'unitless'),
    #### Fuel Grain
    "L_fuel": ToMetric(11.5, 'in'),
    "OD_fuel": ToMetric(3.375, 'in'),
    "Cylindrical": True,
    "ID_fuel": ToMetric(2, 'in'),
    #### Nozzle
    "d_t": ToMetric(1.0, 'in'),
    "alpha": np.deg2rad(15), # Nozzle Diverging Half-Cone Angle
    "d_e": ToMetric(2.22, 'in'),
    ### Ambient Conditions
    "P_amb": ToMetric(12.2, 'psi'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

# Use for Rhode Setup
initialInputs1 = {
    # Purpose:  Dictionary of Initial Inputs, organized by section
    "T_tank": ToMetric(280, 'K'),
    "m_T": ToMetric(2.1, 'kg'),
    "m_N2O": ToMetric(0.3, 'lbm'),
    #### Oxidizer Tank
    "V_tank": ToMetric(0.259, 'L'),
    "P_tank": ToMetric(2900, 'psi'),
    #### Injector
    "A_inj": ToMetric(3.0828e-06, 'm^2'),
    # "A_inj": ToMetric(3.4e-6, 'm^2'),
    "C_d": ToMetric(0.3, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToMetric(2089.83281, 'kg/m^3'),
    "M_chmb": ToMetric(30.169, 'g/mol'),
    "gamma": ToMetric(1.2516, 'unitless'),
    "T_chmb": ToMetric(4000, 'K'),
    ## Fuel Regression Properties
    # "n": ToEnglish(0.3667, 'unitless'),
    # "a": ToEnglish(3.413500729493111e-05, 'unitless'),
    "n": ToMetric(1.6386, 'unitless'),
    "a": ToMetric(9.33E-8, 'unitless'),
    #### Fuel Grain
    "L_fuel": ToMetric(5.0, 'in'),
    "OD_fuel": ToMetric(2, 'in'),
    "Cylindrical": True,
    "ID_fuel": ToMetric(0.5, 'in'),
    #### Nozzle
    "d_t": ToMetric(0.319, 'in'),
    "alpha": np.deg2rad(0), # Nozzle Diverging Half-Cone Angle
    "d_e": ToMetric(0.319, 'in'),
    ### Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

pprint.pprint(initialInputs)

# Creation of the Model, for more information look at Model.py
modelZK = Model(initialInputs, ZK=True)
modelBernoulli = Model(initialInputs, ZK=False)
# Call for the ODE45 Intergration function (solve_ivp)
dfZK = modelZK.ODE45()
dfBe = modelBernoulli.ODE45()

## Formulate CSV files and will put into CSV folder
# print(tabulate(df.iloc[::10], headers='keys', tablefmt='fancy_grid')) # Print every 10th row

totalImpulseZK = dfZK['impulse'].iloc[-1]
totalImpulseBernoulli = dfBe['impulse'].iloc[-1]

totalFinalMassZK = dfZK['mox'].iloc[-1] + dfZK['mf'].iloc[-1]
totalFinalMassBernoulli = dfBe['mox'].iloc[-1] + dfBe['mf'].iloc[-1]

totalMassConsumedZK = ToEnglish(modelZK.mt_i, 'kg') - totalFinalMassZK
totalMassConsumedBernoulli = ToEnglish(modelBernoulli.mt_i, 'kg') - totalFinalMassBernoulli

IspZK = totalImpulseZK/totalMassConsumedZK
IspBernoulli = totalImpulseBernoulli/totalMassConsumedBernoulli

avgThrustZK = totalImpulseZK/dfZK['time'].iloc[-1]
avgThrustBernoulli = totalImpulseBernoulli/dfBe['time'].iloc[-1]

peakThrustZK = dfZK['thrust'].max()
peakThrustBernoulli = dfBe['thrust'].max()

burnTimeZK = dfZK['time'].iloc[-1]
burnTimeBernoulli = dfBe['time'].iloc[-1]

### Console Output ###
print("\n")
print("-------------------- Ullage Fraction -------------------")
print(f'ZK:        {np.round(modelZK.blowdown.ullageFraction, 3)} [unitless]')
print(f'Bernoulli: {np.round(modelBernoulli.blowdown.ullageFraction, 3)} [unitless]')
print("-------------------- Total Mass -------------------")
print(f'ZK:        {np.round(ToEnglish(modelZK.mt_i, "kg"), 3)} [lbm]')
print(f'Bernoulli: {np.round(ToEnglish(modelBernoulli.mt_i, "kg"), 3)} [lbm]')
print("-------------------- Mass Burned -------------------")
print(f'ZK:        {np.round(totalMassConsumedZK, 3)} [lbm]')
print(f'Bernoulli: {np.round(totalMassConsumedBernoulli, 3)} [lbm]')
print("-------------------- Average Thrust -------------------")
print(f'ZK:        {np.round(avgThrustZK, 3)} [lbf]')
print(f'Bernoulli: {np.round(avgThrustBernoulli, 3)} [lbf]')
print("-------------------- Peak Thrust -------------------")
print(f'ZK:        {np.round(peakThrustZK, 3)} [lbf]')
print(f'Bernoulli: {np.round(peakThrustBernoulli, 3)} [lbf]')
print("-------------------- Burn Time -------------------")
print(f'ZK:        {np.round(burnTimeZK, 3)} [s]')
print(f'Bernoulli: {np.round(burnTimeBernoulli, 3)} [s]')
print("-------------------- Total Impulse -------------------")
print(f'ZK:        {np.round(totalImpulseZK, 3)} [lbf-s]')
print(f'Bernoulli: {np.round(totalImpulseBernoulli, 3)} [lbf-s]')
print("-------------------- Specific Impulse -------------------")
print(f'ZK:        {np.round(IspZK, 3)} [s]')
print(f'Bernoulli: {np.round(IspBernoulli, 3)} [s]')
print("\n")

# Oxidizer and Fuel Masses over Time
plot_graph('Rocket Masses Over Time',
            'Time (s)',
            'Mass (lbm)',
            {'y': dfZK['mox'], 'x': dfZK['time'], 'label': 'Oxidizer Mass ZK'},
            {'y': dfBe['mox'], 'x': dfBe['time'], 'label': 'Oxidizer Mass Bernoulli'},
            {'y': dfZK['mf'], 'x': dfZK['time'], 'label': 'Fuel Mass ZK', 'linestyle': '--', 'color': 'r'},
            {'y': dfBe['mf'], 'x': dfBe['time'], 'label': 'Fuel Mass Bernoulli', 'linestyle': ':', 'color': 'g'})

# Rocket Mass Flow Rate over Time
plot_graph('Rocket Mass Flow Rate over Time', 
           'Time (s)', 
           'Mass Flow Rate (lbm/s)',
           {'y': dfZK['dmf'], 'x': dfZK['time'], 'label': 'Fuel Mass Flow ZK'},
           {'y': dfBe['dmf'], 'x': dfBe['time'], 'label': 'Fuel Mass Flow Bernoulli'},
           {'y': dfZK['dmox'], 'x': dfZK['time'], 'label': 'Oxidizer Mass Flow ZK', 'linestyle': '--', 'color': 'r'},
           {'y': dfBe['dmox'], 'x': dfBe['time'], 'label': 'Oxidizer Mass Flow Bernoulli', 'linestyle': ':', 'color': 'g'})

# Impulse vs. Time
plot_graph('Impulse vs. Time', 
           'Time (s)', 
           'Total Impulse Produced (lbf-s)',
           {'y': dfZK['impulse'], 'x': dfZK['time'], 'label': 'ZK'},
           {'y': dfBe['impulse'], 'x': dfBe['time'], 'label': 'Bernoulli'})

# Mixture Ratio vs. Time
plot_graph('Mixture Ratio vs. Time', 
           'Time (s)', 
           'O/F Ratio',
           {'y': dfZK['OF'], 'x': dfZK['time'], 'label': 'ZK'},
           {'y': dfBe['OF'], 'x': dfBe['time'], 'label': 'Bernoulli'})

# Port Radius vs Time
if initialInputs['Cylindrical']:
    noFuelLine = ToEnglish(initialInputs['OD_fuel']/2, 'm')
else:
    noFuelLine = ToEnglish(modelBernoulli.r_final, 'm')
plot_graph('Port Radius vs Time', 
           'Time (s)', 
           'Port Radius (in)',
           {'y': dfZK['r'], 'x': dfZK['time'], 'label': 'Fuel Grain Port Radius ZK'},
           {'y': dfBe['r'], 'x': dfBe['time'], 'label': 'Fuel Grain Port Radius Bernoulli'},
           {'y': noFuelLine, 'label': 'Fuel Grain Outer Diameter', 'linestyle': ':', 'color': 'r'})

# df = pd.read_csv("src/Model/Test_data/01_31_2024_Collected_data/CSV_Files/Burn_2.csv")
# #find index of max thrust
# max_thrust = df['Thrust (lbf)'].max()
# #cut off dataframe from index 0 to max_index
# max_index = df['Thrust (lbf)'].idxmax()
# df = df.iloc[max_index+0:229500]

# # plt.figure()
# # plt.plot(df['Time (s)']-df['Time (s)'].iloc[0], df['N20 Tank Pressure (psi)'], linewidth=0.1, label='Experimental')
# # plt.plot(dfZK['time'], dfZK['Pox'], linewidth=1, label='ZK', linestyle='--')
# # plt.plot(dfBe['time'], dfBe['Pox'], linewidth=1, label='Bernoulli')
# # plt.xlabel('Time (s)')
# # plt.ylabel('Pressure (psi)')
# # plt.title('Tank Pressure vs Time')

# Rocket Pressures Over Time
plot_graph('Rocket Pressures Over Time', 
           'Time (s)', 
           'Pressure (psi)',
           {'y': dfZK['Pc'], 'x': dfZK['time'], 'label': 'Chamber Pressure ZK'},
           {'y': dfBe['Pc'], 'x': dfBe['time'], 'label': 'Chamber Pressure Bernoulli'},
           {'y': dfZK['Pox'], 'x': dfZK['time'], 'label': 'Oxidizer Tank Pressure ZK', 'linestyle': '--', 'color': 'r'},
           {'y': dfBe['Pox'], 'x': dfBe['time'], 'label': 'Oxidizer Tank Pressure Bernoulli', 'linestyle': ':', 'color': 'g'})

# Thrust vs. Time
plot_graph('Thrust vs. Time', 
           'Time (s)', 
           'Thrust (lbf)',
           {'y': dfZK['thrust'], 'x': dfZK['time'], 'label': 'ZK'},
           {'y': dfBe['thrust'], 'x': dfBe['time'], 'label': 'Bernoulli'})
plt.show()

# Define a mapping dictionary from old column names to desired column names
column_mapping = {
    "time": "time (s)",
    "thrust": "thrust (lbf)",
    "thrustN": "thrust (N)",
    "impulse": "impulse (lbf-s)",
    "Pc": "Chamber pressure (psi)",
    "Pox": "Tank pressure (psi)",
    "mox": "Liquid Oxidizer mass (lbm)",
    "mf": "Fuel mass (lbm)",
    "dmox": "Oxidizer mass flow rate (lbm/s)",
    "dmf": "Fuel mass flow rate (lbm/s)",
    "OF": "Oixidizer to Fuel Ratio (unitless)",
    "r": "Port Radius (in)",
    "dm_total_dt": "Total mass flow rate (lbm/s)",
    "cstar": "cstar (in/s)",
    "T_chmb": "Chamber temperature (K)",
    "M_chmb": "Chamber molecular weight (g/mol)",
    "gamma": "gamma (unitless)",
}

# # Format each value in the DataFrame to have 3 decimal places
# dfZK = dfZK.apply(lambda col: col.map(lambda x: format(x, '.3f')))

# # Rename the columns
# dfZK.rename(columns=column_mapping, inplace=True)

# # Save the DataFrame to a CSV file
# dfZK.to_csv("src/Model/CSV/current_data.csv", index=False)