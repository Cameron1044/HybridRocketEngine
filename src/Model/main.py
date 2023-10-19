import numpy as np
import math as m
import matplotlib.pyplot as plt
from tabulate import tabulate
import pandas as pd
import pprint

from Utilities.Model import Model
from Utilities.utilities import ToMetric

initialInputs = { # Dictionary of initial inputs
    # Injector
    "A_inj": ToMetric(9.4E-6, 'm^2'),
    "C_d": ToMetric(0.4, 'unitless'),
    # Fuel Properties
    "rho_fuel": ToMetric(919, 'kg/m^3'),
    "M": ToMetric(23.147, 'g/mol'),
    "T": ToMetric(3015, 'K'),
    "gamma": ToMetric(1.2705, 'unitless'),
    "n": ToMetric(0.681, 'unitless'),
    "a": ToMetric(2.85E-5, 'unitless'),
    # Fuel Grain
    "L_fuel": ToMetric(12, 'in'),
    "OD_fuel": ToMetric(1.688, 'in'),
    "ID_fuel": ToMetric(1.35, 'in'),
    # Oxidizer Tank
    "V_tank": ToMetric(3, 'L'),
    "P_tank": ToMetric(3000, 'psi'),
    "ullage_fraction": ToMetric(0.20, 'unitless'),
    # Nozzle
    "d_t": ToMetric(0.5, 'in'),
    # Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),

    ##### CONSTANTS #####
    "rho_ox": ToMetric(1226, 'kg/m^3'),
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

model = Model(initialInputs)
df = model.ODE45()

print(tabulate(df.iloc[::10], headers='keys', tablefmt='fancy_grid')) # Print every 10th row
df.to_csv("src/Model/CSV/current_data.csv", index=False) # Save dataframe to csv

plt.figure()
plt.grid(True)
plt.plot(df['time'], df['thrust'], linewidth=2)
plt.title('Thrust vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Thrust (lbf)')

plt.figure()
plt.grid(True)
plt.plot(df['time'], df['impulse'], linewidth=2)
plt.title('Impulse vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Total Impulse Produced (lbf-s)')

plt.figure()
plt.grid(True)
plt.plot(df['time'], df['OF'], linewidth=2)
plt.title('Mixture Ratio vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('O/F Ratio')

plt.figure()
plt.plot(df['time'], df['r'], linewidth=2, label="Fuel Grain Port Radius")
plt.axhline(y=initialInputs['OD_fuel']/2/0.0254, color='r', linestyle='--', linewidth=2, label="Fuel Grain Outer Diameter")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend(loc="best")
plt.title("Port Radius vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Port Radius (in)")

plt.figure()
plt.plot(df['time'], df['dmox'], linewidth=2, label="Fuel Mass Flow")
plt.plot(df['time'], df['dmf'], linewidth=2, label="Oxidizer Mass Flow")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend(loc="best")
plt.title("Rocket Mass Flow Rate over Time")
plt.xlabel("Time (s)")
plt.ylabel("Mass Flow Rate (lbm/s)")

plt.figure()
plt.plot(df['time'], df['Pc'], linewidth=2, label="Chamber Pressure")
plt.plot(df['time'], df['Pox'], linewidth=2, label="Oxidizer Tank Pressure")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend(loc="best")
plt.title("Rocket Pressures Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Pressure (psi)")
plt.show()