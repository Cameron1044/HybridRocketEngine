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

# Aluminum model high end prediction
initialInputs = {
    ## Purpose:  Dictionary of Initial Inputs, organized by section
    "ColdFlow": True,
    "T_tank": ToMetric(270, 'K'),
    "m_T": ToMetric(2.1, 'kg'),
    "m_N2O": ToMetric(3, 'lbm'),
    #### Oxidizer Tank
    "V_tank": ToMetric(3, 'L'),
    "P_tank": ToMetric(3000, 'psi'),
    #### Injector
    # "A_inj": ToMetric(1.8e-5, 'm^2'),
    "A_inj": ToMetric(0.038, 'in^2'),
    "C_d": ToMetric(0.5, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToMetric(2089.83281, 'kg/m^3'),
    "M_chmb": ToMetric(30.169, 'g/mol'),
    "gamma": ToMetric(1.2516, 'unitless'),
    "T_chmb": ToMetric(4000, 'K'),
    #### Fuel Regression Properties
    "n": ToMetric(1.681, 'unitless'),
    "a": ToMetric(9.33E-8, 'unitless'),
    #### Fuel Grains
    "L_fuel": ToMetric(11.5, 'in'),
    "OD_fuel": ToMetric(3.375, 'in'),
    "Cylindrical": False,
    ## "ID_fuel": ToMetric(2.75, 'in'),
    #### Nozzle
    "d_t": ToMetric(0.9, 'in'),
    "alpha": np.deg2rad(15), # Nozzle Diverging Half-Cone Angle
    #### Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),
    # "P_amb": ToMetric(827371, 'Pa'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

# Creation of the Model, for more information look at Model.py
modelZK = Model(initialInputs, ZK=True)
modelBernoulli = Model(initialInputs, ZK=False)

# Call for the ODE45 Intergration function (solve_ivp)
dfZK = modelZK.ODE45()
dfBe = modelBernoulli.ODE45()

# initialInputs['C_d'] = ToMetric(0.085, 'unitless')

# Creation of the Model, for more information look at Model.py
modelZK2 = Model(initialInputs, ZK=True)
modelBernoulli2 = Model(initialInputs, ZK=False)

# Call for the ODE45 Intergration function (solve_ivp)
dfZK2 = modelZK2.ODE45()
dfBe2 = modelBernoulli2.ODE45()

# Burn Time
burnTimeZK = dfZK['time'].iloc[-1]
burnTimeBernoulli = dfBe['time'].iloc[-1]

### Console Output ###
print("-------------------- Ullage Fraction -------------------")
print(f'ZK:        {np.round(modelZK.blowdown.ullageFraction, 3)} [unitless]')
print(f'Bernoulli: {np.round(modelBernoulli.blowdown.ullageFraction, 3)} [unitless]')
print("-------------------- Burn Time -------------------")
print(f'ZK:        {np.round(burnTimeZK, 3)} [s]')
print(f'Bernoulli: {np.round(burnTimeBernoulli, 3)} [s]')
print("\n")

df1 = pd.read_csv("src/Model/Test_data/rhode/cold1.csv")
df2 = pd.read_csv("src/Model/Test_data/rhode/cold2.csv")
df3 = pd.read_csv("src/Model/Test_data/rhode/cold3.csv")
# max_thrust = df1['N20 Tank PT (V)'].max()
# max_index = df1['N20 Tank PT (V)'].idxmax()
# SR = 1/1700
# offset1 = int(115.4/SR)
# offset2 = int(30/SR) + offset1
# df = df1.iloc[max_index + offset1:max_index + offset2]
# offset1 = int((115.3+4.85)/SR)
# offset2 = int(30/SR) + offset1
# df = df1.iloc[max_index + offset1:max_index + offset2]
# df = df1.iloc[0:-1]

plt.figure()
plt.plot(df1['Time [s]']-df1['Time [s]'].iloc[0], df1['Pressure0 [psi]'], linewidth=0.5, label='Coldflow 1')
plt.plot(df2['Time [s]']-df2['Time [s]'].iloc[0], df2['Pressure0 [psi]'], linewidth=0.5, label='Coldflow 2')
plt.plot(df3['Time [s]']-df3['Time [s]'].iloc[0], df3['Pressure0 [psi]'], linewidth=0.5, label='Coldflow 3')
# plt.axvline(x=2.5, linewidth=1, label="Point of Blowout", linestyle='--', color='r')
plt.legend(loc="best")
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.title('Tank Pressure vs Time')


plt.figure()
plt.plot(df1['Time [s]']-df1['Time [s]'].iloc[0], df1['Pressure0 [psi]'], linewidth=0.5, label='Coldflow 1')
plt.plot(df2['Time [s]']-df2['Time [s]'].iloc[0], df2['Pressure0 [psi]'], linewidth=0.5, label='Coldflow 2')
plt.plot(df3['Time [s]']-df3['Time [s]'].iloc[0], df3['Pressure0 [psi]'], linewidth=0.5, label='Coldflow 3')
plt.plot(dfBe['time'], dfBe['Pox'], linewidth=2, label='Bernoulli Cd=0.035', color='orange')
plt.plot(dfZK['time'], dfZK['Pox'], linewidth=1, label='ZK Cd=0.035', linestyle='--', color='red')
# plt.plot(dfBe2['time'], dfBe2['Pox'], linewidth=2, label='Bernoulli Cd=0.085', color='green')
# plt.plot(dfZK2['time'], dfZK2['Pox'], linewidth=1, label='ZK Cd=0.085', linestyle='--', color='purple')
plt.legend(loc="best")
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.title('Tank Pressure vs Time')
plt.show()
# plt.figure()
# plt.plot(df['Time']-df['Time'].iloc[0], df['Combustion Chamber PT (V)']*5000/4.5, linewidth=0.5, label='Experimental')
# plt.plot(dfBe['time'], dfBe['Pc'], linewidth=2, label='Bernoulli Cd=0.035', color='orange')
# plt.plot(dfZK['time'], dfZK['Pc'], linewidth=1, label='ZK Cd=0.035', linestyle='--', color='red')
# # plt.axvline(x=2.5, linewidth=1, label="Point of Blowout", linestyle='--', color='r')
# plt.legend(loc="best")
# plt.xlabel('Time (s)')
# plt.ylabel('Pressure (psi)')
# plt.title('Tank Pressure vs Time')

# plt.figure()
# plt.plot(df['Time']-df['Time'].iloc[0], df['Thrust LC1 (lbf)']*-1 + df['Thrust LC3 (lbf)'], linewidth=0.5, label='Experimental')
# plt.plot(dfBe['time'], dfBe['thrust'], linewidth=2, label='Bernoulli', color='orange')
# plt.plot(dfZK['time'], dfZK['thrust'], linewidth=1, label='ZK', linestyle='--', color='g')
# plt.legend(loc="best")
# plt.xlabel('Time (s)')
# plt.ylabel('Thrust (lbf)')
# plt.title('Thrust vs Time')

# plt.figure()
# plt.plot(df1['Time']-df1['Time'].iloc[0], df1['Thrust LC1 (lbf)']*-1 + df1['Thrust LC3 (lbf)'], linewidth=0.5, label='Experimental')
# plt.legend(loc="best")
# plt.xlabel('Time (s)')
# plt.ylabel('Thrust (lbf)')
# plt.title('Thrust vs Time')

# plt.figure()
# plt.plot(df['Time']-df['Time'].iloc[0], (df['Combustion Chamber PT (V)'] - 0.498)*5000/4.5, linewidth=0.5, label='Experimental')
# plt.legend(loc="best")
# plt.xlabel('Time (s)')
# plt.ylabel('Combustion Chamber (psi)')
# plt.title('Thrust vs Time')

# plt.figure()
# ax1 = plt.gca()
# ax1.plot(df1['Time'] - df1['Time'].iloc[0], (df1['N20 Tank PT (V)'])*5000/4.5, linewidth=0.5, label='Ox Tank Pressure', color='b')
# ax1.set_xlabel('Time (s)')
# ax1.set_ylabel('N20 Tank PT (V)', color='b')
# ax1.tick_params(axis='y', labelcolor='b')
# ax2 = ax1.twinx()
# ax2.plot(df1['Time'] - df1['Time'].iloc[0], df1['Thrust LC1 (lbf)']*-1 + df1['Thrust LC2 (lbf)'] + df1['Thrust LC3 (lbf)'] + 85, linewidth=0.5, label='Thrust', color='r')
# ax2.set_ylabel('Thrust (lbf)', color='r')
# ax2.tick_params(axis='y', labelcolor='r')
# ax1.legend(loc='upper left')
# ax2.legend(loc='upper right')
# plt.title('N20 Tank PT and Thrust over Time')
# plt.show()

# # plt.figure()
# # plt.plot(df1['Fluids PT1(psi)'], linewidth=1, label='Experimental')
# # plt.plot(df1['Fluids PT2 (psi)'], linewidth=1, label='Experimental')
# plt.figure()
# plt.plot(df['Time']-df['Time'].iloc[0], df['Combustion Chamber PT (V)']*5000/4.5, linewidth=1, label='Experimental')
# plt.show()

# plt.plot(df1['Time']-df1['Time'].iloc[0], df1['N20 Tank PT (V)']*5000/4.5, linewidth=0.5, label='Experimental')