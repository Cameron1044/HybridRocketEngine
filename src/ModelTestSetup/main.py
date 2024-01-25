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

# initialInputs = {
#     # Purpose:  Dictionary of Initial Inputs, organized by section
#     "T_tank": ToMetric(288, 'K'),
#     "m_T": ToMetric(2.1, 'kg'),
#     "m_N2O": ToMetric(2.1, 'kg'),
#     #### Oxidizer Tank
#     "V_tank": ToMetric(3, 'L'),
#     "P_tank": ToMetric(3000, 'psi'),
#     #### Injector
#     "A_inj": ToMetric(1.8e-5, 'm^2'),
#     "C_d": ToMetric(0.4, 'unitless'),
#     #### Fuel Properties
#     "rho_fuel": ToMetric(2089.83281, 'kg/m^3'),
#     ## Fuel Regression Properties
#     "n": ToMetric(1.681, 'unitless'),
#     "a": ToMetric(9.33E-8, 'unitless'),
#     #### Fuel Grain
#     "L_fuel": ToMetric(12, 'in'),
#     "OD_fuel": ToMetric(3.375, 'in'),
#     #### Nozzle
#     "d_t": ToMetric(1.2, 'in'),
#     "alpha": np.deg2rad(15), # Nozzle Diverging Half-Cone Angle
#     ### Ambient Conditions
#     "P_amb": ToMetric(102675.3, 'Pa'),
#     ##### CONSTANTS #####
#     "Ru": ToMetric(8.3143, 'J/(mol*K)'),
# }
initialInputs = {
    # Purpose:  Dictionary of Initial Inputs, organized by section
    "T_tank": ToMetric(288, 'K'),
    "m_T": ToMetric(2.1, 'kg'),
    "m_N2O": ToMetric(2.1, 'kg'),
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
    "a": ToEnglish(0.104, 'unitless'),
    #### Fuel Grain
    "L_fuel": ToMetric(11.5, 'in'),
    "OD_fuel": ToMetric(3.375, 'in'),
    "ID_fuel": ToMetric(2, 'in'),
    #### Nozzle
    "d_t": ToMetric(1.0, 'in'),
    "alpha": np.deg2rad(15), # Nozzle Diverging Half-Cone Angle
    ### Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

initialInputs2 = {
    # Purpose:  Dictionary of Initial Inputs, organized by section
    "T_tank": ToMetric(288, 'K'),
    "m_T": ToMetric(2.1, 'kg'),
    "m_N2O": ToMetric(0.185, 'kg'),
    #### Oxidizer Tank
    "V_tank": ToMetric(0.259, 'L'),
    "P_tank": ToMetric(2300, 'psi'),
    #### Injector
    # "A_inj": ToMetric(1.8e-6, 'm^2'),
    "A_inj": ToMetric(3.4e-6, 'm^2'),
    "C_d": ToMetric(0.175, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToMetric(2089.83281, 'kg/m^3'),
    ## Fuel Regression Properties
    "n": ToMetric(1.681, 'unitless'),
    "a": ToMetric(9.33E-8, 'unitless'),
    #### Fuel Grain
    "L_fuel": ToMetric(5.5, 'in'),
    "OD_fuel": ToMetric(2, 'in'),
    "ID_fuel": ToMetric(0.5, 'in'),
    #### Nozzle
    "d_t": ToMetric(0.346, 'in'),
    "alpha": np.deg2rad(15), # Nozzle Diverging Half-Cone Angle
    ### Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

initialInputs1 = {
    # Purpose:  Dictionary of Initial Inputs, organized by section
    "T_tank": ToEnglish(288, 'K'),
    "m_T": ToEnglish(2.1, 'kg'),
    "m_N2O": ToEnglish(0.185, 'kg'),
    #### Oxidizer Tank
    "V_tank": ToEnglish(0.259, 'L'),
    "P_tank": ToEnglish(2300, 'psi'),
    #### Injector
    "A_inj": ToEnglish(1.9e-6, 'm^2'),
    "C_d": ToEnglish(0.175, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToEnglish(2089.83281, 'kg/m^3'),
    "M_chmb": ToEnglish(23.719/1000, 'unitless'),
    "gamma": ToEnglish(1.2641, 'unitless'),
    "T_chmb": ToEnglish(3075, 'K'),
    ## Fuel Regression Properties
    "n": ToEnglish(1.681, 'unitless'),
    "a": ToEnglish(9.33E-8/(0.0254**(1 + 2*(1.681)))*(0.453592**(-1.681)), 'unitless'),
    #### Fuel Grain
    "L_fuel": ToEnglish(5.5, 'in'),
    "OD_fuel": ToEnglish(2, 'in'),
    "ID_fuel": ToEnglish(0.5, 'in'),
    #### Nozzle
    "d_t": ToEnglish(0.4, 'in'),
    "alpha": np.deg2rad(27), # Nozzle Diverging Half-Cone Angle
    ### Ambient Conditions
    "P_amb": ToEnglish(102675.3, 'Pa'),
    ##### CONSTANTS #####
    # "Ru": ToEnglish(8.3143, 'J/(mol*K)'),
}

pprint.pprint(initialInputs)

# "rho_fuel": ToMetric(919, 'kg/m^3'),
# ## Fuel Regression Properties
# "n": ToMetric(0.681, 'unitless'),
# "a": ToMetric(2.85E-5, 'unitless'),
# initialInputs['T'] = ToMetric(4000, 'K')
# initialInputs['M'] = ToMetric(30.169, 'g/mol')
# initialInputs['gamma'] = ToMetric(1.2516, 'unitless')
# initialInputs['rho_fuel'] = ToMetric(1994, 'kg/m^3')
# initialInputs['n'] = ToMetric(1.681, 'unitless')
# initialInputs['a'] = ToMetric(9.33E-8, 'unitless')
# initialInputs['A_inj'] = ToMetric(6.4e-6, 'm^2')
# initialInputs['d_t'] = ToMetric(0.5, 'in')

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

### PLOTTING ####
# plot_graph('Combustion Chamber Temperature Over Time',
#             'Time (s)',
#             'Temperature (K)',
#             {'y': dfZK['T_chmb'], 'x': dfZK['time'], 'label': 'Chamber Temperature ZK'},
#             {'y': dfBe['T_chmb'], 'x': dfBe['time'], 'label': 'Chamber Temperature Bernoulli'})

# plot_graph('Combustion Chamber Molecular Weight Over Time',
#             'Time (s)',
#             'Molecular Weight (g/mol)',
#             {'y': dfZK['M_chmb'], 'x': dfZK['time'], 'label': 'Chamber Molecular Weight ZK'},
#             {'y': dfBe['M_chmb'], 'x': dfBe['time'], 'label': 'Chamber Molecular Weight Bernoulli'})

# plot_graph('Combustion Chamber Gamma Over Time',
#             'Time (s)',
#             'Gamma (unitless)',
#             {'y': dfZK['gamma'], 'x': dfZK['time'], 'label': 'Gamma ZK'},
#             {'y': dfBe['gamma'], 'x': dfBe['time'], 'label': 'Gamma Bernoulli'})

# # C* vs. Time
# plot_graph('C* vs. Time', 
#            'Time (s)', 
#            'C* (in/s)',
#            {'y': dfZK['cstar'], 'x': dfZK['time'], 'label': 'ZK'},
#            {'y': dfBe['cstar'], 'x': dfBe['time'], 'label': 'Bernoulli'})

# # Oxidizer and Fuel Masses over Time
# plot_graph('Rocket Masses Over Time',
#             'Time (s)',
#             'Mass (lbm)',
#             {'y': dfZK['mox'], 'x': dfZK['time'], 'label': 'Oxidizer Mass ZK'},
#             {'y': dfBe['mox'], 'x': dfBe['time'], 'label': 'Oxidizer Mass Bernoulli'},
#             {'y': dfZK['mf'], 'x': dfZK['time'], 'label': 'Fuel Mass ZK', 'linestyle': '--', 'color': 'r'},
#             {'y': dfBe['mf'], 'x': dfBe['time'], 'label': 'Fuel Mass Bernoulli', 'linestyle': ':', 'color': 'g'})

# # Rocket Mass Flow Rate over Time
# plot_graph('Rocket Mass Flow Rate over Time', 
#            'Time (s)', 
#            'Mass Flow Rate (lbm/s)',
#            {'y': dfZK['dmf'], 'x': dfZK['time'], 'label': 'Fuel Mass Flow ZK'},
#            {'y': dfBe['dmf'], 'x': dfBe['time'], 'label': 'Fuel Mass Flow Bernoulli'},
#            {'y': dfZK['dmox'], 'x': dfZK['time'], 'label': 'Oxidizer Mass Flow ZK', 'linestyle': '--', 'color': 'r'},
#            {'y': dfBe['dmox'], 'x': dfBe['time'], 'label': 'Oxidizer Mass Flow Bernoulli', 'linestyle': ':', 'color': 'g'})

# Mixture Ratio vs. Time
plot_graph('Mixture Ratio vs. Time', 
           'Time (s)', 
           'O/F Ratio',
           {'y': dfZK['OF'], 'x': dfZK['time'], 'label': 'ZK'},
           {'y': dfBe['OF'], 'x': dfBe['time'], 'label': 'Bernoulli'})

# Port Radius vs Time
plot_graph('Port Radius vs Time', 
           'Time (s)', 
           'Port Radius (in)',
           {'y': dfZK['r'], 'x': dfZK['time'], 'label': 'Fuel Grain Port Radius ZK'},
           {'y': dfBe['r'], 'x': dfBe['time'], 'label': 'Fuel Grain Port Radius Bernoulli'},
           {'y': ToEnglish(initialInputs['OD_fuel']/2, 'm'), 'label': 'Fuel Grain Outer Diameter', 'linestyle': ':', 'color': 'r'})

# # Impulse vs. Time
# plot_graph('Impulse vs. Time', 
#            'Time (s)', 
#            'Total Impulse Produced (lbf-s)',
#            {'y': dfZK['impulse'], 'x': dfZK['time'], 'label': 'ZK'},
#            {'y': dfBe['impulse'], 'x': dfBe['time'], 'label': 'Bernoulli'})

# Thrust vs. Time
plot_graph('Thrust vs. Time', 
           'Time (s)', 
           'Thrust (lbf)',
           {'y': dfZK['thrust'], 'x': dfZK['time'], 'label': 'ZK'},
           {'y': dfBe['thrust'], 'x': dfBe['time'], 'label': 'Bernoulli'})

# Rocket Pressures Over Time
plot_graph('Rocket Pressures Over Time', 
           'Time (s)', 
           'Pressure (psi)',
           {'y': dfZK['Pc'], 'x': dfZK['time'], 'label': 'Chamber Pressure ZK'},
           {'y': dfBe['Pc'], 'x': dfBe['time'], 'label': 'Chamber Pressure Bernoulli'},
           {'y': dfZK['Pox'], 'x': dfZK['time'], 'label': 'Oxidizer Tank Pressure ZK', 'linestyle': '--', 'color': 'r'},
           {'y': dfBe['Pox'], 'x': dfBe['time'], 'label': 'Oxidizer Tank Pressure Bernoulli', 'linestyle': ':', 'color': 'g'})

# Creation of the Model, for more information look at Model.py
# model1 = Model(initialInputs, ZK=True)
# initialInputs['C_d'] = ToMetric(0.2, 'unitless')
# model2 = Model(initialInputs, ZK=True)
# initialInputs['C_d'] = ToMetric(0.6, 'unitless')
# model3 = Model(initialInputs, ZK=True)
# # Call for the ODE45 Intergration function (solve_ivp)
# df1 = model1.ODE45()
# df2 = model2.ODE45()
# df3 = model3.ODE45()

# plot_graph('Thrust vs. Time', 
#            'Time (s)', 
#            'Pressure (psi)',
#            {'y': df1['thrust'], 'x': df1['time'], 'label': 'Cd = 0.4', 'color': 'r'},
#            {'y': df2['thrust'], 'x': df2['time'], 'label': 'Cd = 0.2', 'color': 'g'},
#            {'y': df3['thrust'], 'x': df3['time'], 'label': 'Cd = 0.6', 'color': 'b'})

# plot_graph('Pressures vs. Time', 
#            'Time (s)', 
#            'Pressure (psi)',
#            {'y': df1['Pc'], 'x': df1['time'], 'label': 'Pc: Cd = 0.4', 'linestyle': '--', 'color': 'r'},
#            {'y': df2['Pc'], 'x': df2['time'], 'label': 'Pc: Cd = 0.2', 'linestyle': '--', 'color': 'g'},
#            {'y': df3['Pc'], 'x': df3['time'], 'label': 'Pc: Cd = 0.6', 'linestyle': '--', 'color': 'b'},
#            {'y': df1['Pox'], 'x': df1['time'], 'label': 'Pox: Cd = 0.4', 'color': 'r'},
#            {'y': df2['Pox'], 'x': df2['time'], 'label': 'Pox: Cd = 0.2', 'color': 'g'},
#            {'y': df3['Pox'], 'x': df3['time'], 'label': 'Pox: Cd = 0.6', 'color': 'b'},)

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

# Format each value in the DataFrame to have 3 decimal places
dfZK = dfZK.apply(lambda col: col.map(lambda x: format(x, '.3f')))

# Rename the columns
dfZK.rename(columns=column_mapping, inplace=True)

# Save the DataFrame to a CSV file
dfZK.to_csv("src/Model/CSV/current_data.csv", index=False)