import numpy as np
import math as m
import matplotlib.pyplot as plt
from tabulate import tabulate
import pandas as pd
import pprint

from Utilities.Model import Model
from Utilities.utilities import ToMetric, ToEnglish

initialInputs = {
    # Purpose:  Dictionary of Initial Inputs, organized by section
    
    #### Oxidizer Tank
    "V_tank": ToMetric(3, 'L'),
    "P_tank": ToMetric(3000, 'psi'),
    "ullage_fraction": ToMetric(0.2, 'unitless'),
    "rho_ox": ToMetric(1226, 'kg/m^3'),
    #### Injector
    "A_inj": ToMetric(9.4e-6, 'm^2'),
    "C_d": ToMetric(0.4, 'unitless'),
    #### Fuel Properties
    "rho_fuel": ToMetric(919, 'kg/m^3'),
    ## Fuel Regression Properties
    "n": ToMetric(0.681, 'unitless'),
    "a": ToMetric(2.85E-5, 'unitless'),
    ## Combustion Fuel Properties
    # "gamma": ToMetric(1.2705, 'unitless'),
    # "M": ToMetric(23.147, 'g/mol'),
    # "T": ToMetric(3015, 'K'),
    "gamma": ToMetric(1.2448, 'unitless'),
    "M": ToMetric(27.442, 'g/mol'),
    "T": ToMetric(3474, 'K'),
    #### Fuel Grain
    "L_fuel": ToMetric(6, 'in'),
    "OD_fuel": ToMetric(3.375, 'in'),
    "ID_fuel": ToMetric(1.0, 'in'),
    #### Nozzle
    "d_t": ToMetric(0.5, 'in'),
    #### Ambient Conditions
    "P_amb": ToMetric(102675.3, 'Pa'),
    ##### CONSTANTS #####
    "Ru": ToMetric(8.3143, 'J/(mol*K)'),
}

# initialInputs['T'] = ToMetric(4000, 'K')
# initialInputs['M'] = ToMetric(30.169, 'g/mol')
# initialInputs['gamma'] = ToMetric(1.2516, 'unitless')
# initialInputs['rho_fuel'] = ToMetric(1994, 'kg/m^3')
# initialInputs['n'] = ToMetric(1.681, 'unitless')
# initialInputs['a'] = ToMetric(9.33E-8, 'unitless')
# initialInputs['A_inj'] = ToMetric(6.4e-6, 'm^2')
# initialInputs['d_t'] = ToMetric(0.5, 'in')

# Creation of the Model, for more information look at Model.py
model = Model(initialInputs)
# Call for the ODE45 Intergration function (solve_ivp)
df = model.ODE45()

#### PLOTTING ####
## Formulate CSV files and will put into CSV folder
# print(tabulate(df.iloc[::10], headers='keys', tablefmt='fancy_grid')) # Print every 10th row
df.to_csv("src/Model/CSV/current_data.csv", index=False) # Save dataframe to csv

moi = initialInputs['rho_ox']*initialInputs['V_tank']*0.2
total_impulse = df['impulse'].iloc[-1]

Isp = total_impulse/ToEnglish(model.mt_i, 'kg')
# print(ToEnglish(model.mt_i, 'kg'))
print(f'Total Impulse: {total_impulse} [lbm-s]')
print(f'Specific Impulse: {Isp} [s]')

# # ## Plotting Values from df model
# # # Thrust Over Time
# plt.figure()
# plt.grid(True)
# plt.plot(df['time'], df['cstar'], linewidth=2)
# plt.title('Thrust vs. Time')
# plt.xlabel('Time (s)')
# plt.ylabel('Thrust (lbf)')

## Plotting Values from df model
# Thrust Over Time
plt.figure()
plt.plot(df['time'], df['thrust'], linewidth=2)
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.title('Thrust vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Thrust (lbf)')

# Impulse Over Time
plt.figure()
plt.plot(df['time'], df['impulse'], linewidth=2)
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.title('Impulse vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Total Impulse Produced (lbf-s)')

# OF Ratio Over Time
plt.figure()
plt.plot(df['time'], df['OF'], linewidth=2)
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.title('Mixture Ratio vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('O/F Ratio')

# Fuel Grain Port Radius Over Time
plt.figure()
plt.plot(df['time'], df['r'], linewidth=2, label="Fuel Grain Port Radius")
plt.axhline(y=initialInputs['OD_fuel']/2/0.0254, color='r', linestyle='--', linewidth=2, label="Fuel Grain Outer Diameter")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend(loc="best")
plt.title("Port Radius vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Port Radius (in)")

# Rocket Mass Flow Rates Over Time
plt.figure()
plt.plot(df['time'], df['dmf'], linewidth=2, label="Fuel Mass Flow")
plt.plot(df['time'], df['dmox'], linewidth=2, label="Oxidizer Mass Flow")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend(loc="best")
plt.title("Rocket Mass Flow Rate over Time")
plt.xlabel("Time (s)")
plt.ylabel("Mass Flow Rate (lbm/s)")

# Rocket Pressure Over Time
plt.figure()
plt.plot(df['time'], df['Pc'], linewidth=2, label="Chamber Pressure")
plt.plot(df['time'], df['Pox'], linewidth=2, label="Oxidizer Tank Pressure")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend(loc="best")
plt.title("Rocket Pressures Over Time")
plt.xlabel("Time (s)")
plt.ylabel("Pressure (psi)")

# Define a range of injector areas for analysis
a_inj_range = np.linspace(5E-6, 2E-5, 8)

# Plot thrust vs time for different injector areas
plt.figure()
for a_inj in a_inj_range:
    # Update the injector area in the initial inputs
    initialInputs['A_inj'] = ToMetric(a_inj, 'm^2')
    
    # Create a model instance and solve the ODEs
    model = Model(initialInputs)
    df = model.ODE45()
    
    # Plot the results
    plt.plot(df['time'], df['thrust'], label=f"{a_inj:.2E} m^2", linewidth=1)
    plt.legend(loc="best")
    plt.grid(which='both', linestyle='--', linewidth=0.5)
    plt.title("Thrust vs Time for Various Injector Areas")
    plt.xlabel("Time (s)")
    plt.ylabel("Thrust (lbf)")

# Plot chamber pressure vs time for different injector areas
plt.figure()
for a_inj in a_inj_range:
    initialInputs['A_inj'] = ToMetric(a_inj, 'm^2')
    model = Model(initialInputs)
    df = model.ODE45()
    plt.plot(df['time'], df['Pc'], label=f"{a_inj:.2E} m^2", linewidth=1)
    plt.legend(loc="best")
    plt.grid(which='both', linestyle='--', linewidth=0.5)
    plt.title("Chamber Pressure vs Time for Various Injector Areas")
    plt.xlabel("Time (s)")
    plt.ylabel("Chamber Pressure (psi)")

# Reset the injector area to its original value
initialInputs['A_inj'] = ToMetric(9.4e-6, 'm^2')

# Define a range of throat diameters for analysis
d_t_range = np.linspace(0.4, 0.7, 8)

# Plot thrust vs time for different throat diameters
plt.figure()
for d_t in d_t_range:
    # Update the throat diameter in the initial inputs
    initialInputs['d_t'] = ToMetric(d_t, 'in')
    
    # Create a model instance and solve the ODEs
    model = Model(initialInputs)
    df = model.ODE45()
    
    # Plot the results
    plt.plot(df['time'], df['thrust'], label=f"{d_t:.2} in", linewidth=1)
    plt.legend(loc="best")
    plt.grid(which='both', linestyle='--', linewidth=0.5)
    plt.title("Thrust vs Time for Various Throat Diamaters")
    plt.xlabel("Time (s)")
    plt.ylabel("Thrust (lbf)")

# Plot chamber pressure vs time for different throat diameters
plt.figure()
for d_t in d_t_range:
    initialInputs['d_t'] = ToMetric(d_t, 'in')
    model = Model(initialInputs)
    df = model.ODE45()
    plt.plot(df['time'], df['Pc'], label=f"{d_t:.2} in", linewidth=1)
    plt.legend(loc="best")
    plt.grid(which='both', linestyle='--', linewidth=0.5)
    plt.title("Chamber Pressure vs Time for Various Throat Diamaters")
    plt.xlabel("Time (s)")
    plt.ylabel("Chamber Pressure (psi)")

# Reset the throat diameter to its original value
initialInputs['d_t'] = ToMetric(0.5, 'in')

# Create model instances for different fuel mixtures
model1 = Model(initialInputs)
# initialInputs['T'] = ToMetric(5583.19, 'K')
# initialInputs['M'] = ToMetric(39.992, 'g/mol')
# initialInputs['gamma'] = ToMetric(1.1587, 'unitless')
# initialInputs['rho_fuel'] = ToMetric(2767.99, 'kg/m^3')
initialInputs['T'] = ToMetric(4000, 'K')
initialInputs['M'] = ToMetric(30.169, 'g/mol')
initialInputs['gamma'] = ToMetric(1.2516, 'unitless')
initialInputs['rho_fuel'] = ToMetric(1994, 'kg/m^3')
initialInputs['n'] = ToMetric(1.681, 'unitless')
initialInputs['a'] = ToMetric(9.33E-8, 'unitless')
# initialInputs['A_inj'] = ToMetric(6.4e-6, 'm^2')
# initialInputs['d_t'] = ToMetric(0.5, 'in')
model2 = Model(initialInputs)

# Solve the ODEs for each model
df1 = model1.ODE45()
df2 = model2.ODE45()

# Plot thrust vs time for different fuel mixtures
plt.figure()
plt.grid(True)
plt.plot(df1['time'], df1['thrust'], linewidth=2, label="0% Aluminum")
plt.plot(df2['time'], df2['thrust'], linewidth=2, label="60% Aluminum")
plt.title('Thrust vs. Time for Various Fuel Mixtures')
plt.xlabel('Time (s)')
plt.ylabel('Thrust (lbf)')

# Plot chamber pressure vs time for different fuel mixtures
plt.figure()
plt.plot(df1['time'], df1['Pc'], linewidth=2, label="0% Aluminum")
plt.plot(df2['time'], df2['Pc'], linewidth=2, label="60% Aluminum")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend(loc="best")
plt.title("Chamber Pressures Over Time for Various Fuel Mixtures")
plt.xlabel("Time (s)")
plt.ylabel("Pressure (psi)")

# Plot tank pressure vs time for different fuel mixtures
plt.figure()
plt.plot(df1['time'], df1['Pox'], linewidth=2, label="0% Aluminum")
plt.plot(df2['time'], df2['Pox'], linewidth=2, label="100% Aluminum")
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.legend(loc="best")
plt.title("Tank Pressures Over Time for Various Fuel Mixtures")
plt.xlabel("Time (s)")
plt.ylabel("Pressure (psi)")

# Display all the plots
plt.show()
