import numpy as np
import math as m
import matplotlib.pyplot as plt
from tabulate import tabulate
import pandas as pd

from Utilities.blowdownModel import BlowdownModel
from Utilities.utilities import units, plot_as_individual, plot_as_subplots

def mainLoop(tf, tstep, initial_inputs):
    "----- Inital Conditions -----"
    # Injector
    Ainj = initial_inputs["Ainj"]
    Cd = initial_inputs["Cd"]
    # Fuel Properties
    rho_fuel = initial_inputs["rho_fuel"]
    M = initial_inputs["M"]
    T = initial_inputs["T"]
    gamma = initial_inputs["gamma"]
    n = initial_inputs["n"]
    a = initial_inputs["a"]
    # Fuel Grain
    L = initial_inputs["L"]
    r_port0 = initial_inputs["r_port0"]
    # Oxidizer Tank
    m_loaded = initial_inputs["m_loaded"]
    n_gas = initial_inputs["n_gas"]
    V_tank = initial_inputs["V_tank"]
    To = initial_inputs["Ti_tank"]
    m_T = initial_inputs["m_T"]
    # Nozzle
    d_t = initial_inputs["d_t"]
    cd_throat = initial_inputs["cd_throat"]
    # Volumes
    V_chmb_emty = initial_inputs["V_chmb_emty"]

    "----- Initial Calculations -----"
    # INITIAL COMBUSTION CHAMBER
    A_t = np.pi*(d_t/2)**2                                          # Cross-sectional area of nozzle throat | [m^2]
    A_port = np.pi * r_port0**2                                # Inital fuel grain burn surface area | [m^2]
    r = r_port0                                                     # Initial fuel grain radius | [m]

    "----- Initial Setup -----"
    chem = BlowdownModel(gas="Ar", mass_loaded=m_loaded, tank_volume=V_tank, initial_temp=Ti_tank)
    n_go = chem.n_go
    n_lo = chem.n_lo
    P_chmb = 0

    # Initial derivatives
    dr_dt = 0
    dm_ox_dt = 0
    dm_f_dt = 0

    Fuel_OD = units(1.688, "in_to_m")       # Fuel Grain Outer Diameter | [in --> m] 
    Fuel_ID = units(1.25, "in_to_m")        # Fuel Grain Inner Diameter | [in --> m]
    m_ox = m_loaded
    m_fuel = rho_fuel * np.pi * (Fuel_OD**2-Fuel_ID**2)/4

    # time Loop
    i_f=tf/tstep

    for i in range(0,int(i_f)):
        "----- BLOWDOWN -----"
        dT, dn_g, dn_l, P = chem.ZKModel(To, n_go, n_lo, n_gas, V_tank, Ainj, Cd, P_chmb, m_T)
        # Forward Difference Method
        To = To + dT*tstep
        n_go = n_go + dn_g*tstep
        n_lo = n_lo + dn_l*tstep
        m_ox = m_ox + dm_ox_dt*tstep                                               # [kg]
        m_fuel = m_fuel + dm_f_dt*tstep                                            # [kg]

        "----- COMBUSTION CHAMBER -----"
        # Calculating Mass flow of Liquid Oxidizer
        dm_ox_dt = -dn_l*chem.MW_N2O                                      # [kg/s]

        # Calculating change in fuel regression radius
        dr_dt = -a * (dm_ox_dt / A_port)**n                                # [m/s]
        r = r + dr_dt*tstep                                                 # [m]

        A_port = np.pi * r**2                                        # [m^2]
        # Calculating Fuel Burn Surface Area
            # Note: Cylindrical Surface area --> 2 pi r h
        A_port_surf = 2 * np.pi * r * L                                          # [m^2]

        # Calculating mass of fuel flow rate
        dm_f_dt = -dr_dt * A_port_surf * rho_fuel                                 # [kg/s] 

        # Calculating Oxidizer to Fuel Ratio
        OF = dm_ox_dt / dm_f_dt

        # Calculating Chamber Pressure
        C_star_nom = np.sqrt((chem.R * T) / (gamma * M))                         # [m/s]
        C_star_denom = (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))

            # Note: This relation shows up in mdot exit and C_star | Dimensionless
        C_star = C_star_nom / C_star_denom  

        # Calculating Chamber Pressure
        P_chmb = ((dm_ox_dt + dm_f_dt) * C_star)/(A_t)                      # [Pa]

        # Calculating Chamber mass flow rate
        m_chmb = m_fuel + m_ox                                              # [kg/s]

        # Calculating new volume in the Combusation Chamber
        V_chmb = V_chmb_emty - (np.pi * r**2 * L)                           # [m^3]

        # Calculating Density of the combustion chamber
        rho_chmb = m_chmb / V_chmb                                          # [kg/m^3]

        "----- NOZZLE -----"
        # Calculating the choked mass flow of the combustion chamber
        CompressibleFactor = (2 / (gamma + 1))**((gamma + 1) / (gamma - 1))
        mdot_exit = cd_throat * A_t * np.sqrt(gamma * rho_chmb * P_chmb * CompressibleFactor) # [kg/s]
        # Calculating exit velocity
        v_exit = np.sqrt(((2*gamma)/(gamma-1)) * ((chem.R*T)/M))
            # Note: Assumption for perfectly expanded nozzle such that P Chamber = P Exit
        # Calculating thrust
        thrust = mdot_exit * v_exit                                         # [N]
        # If chamber exceeds oxidizer tank pressure stop loop.

        "----- DATAFRAME -----"
        new_data = { # Dictionary of new data
            'time': i * tstep,
            'dm_ox_dt': dm_ox_dt,
            'dr_dt': dr_dt,
            'thrust': units(thrust, "N_to_lbf"),
            'n_lo': n_lo,
            'n_go': n_go,
            'OF': OF,
            'r': units(r, "m_to_in"),
            'V_chmb': units(V_chmb, "m3_to_in3"),
            'p_chmb': units(P_chmb, "pa_to_psi"),
            'p_tank': units(P, "pa_to_psi")
        }
        df.loc[len(df)] = new_data # Append new data to dataframe

        "----- BREAK CONDITIONS -----"
        tf1 = P_chmb > P # Chamber pressure is greater than tank pressure
        tf2 = n_lo <= 0 # No more liquid oxidizer
        tf3 = np.isnan(P) or np.isnan(P_chmb) or np.isnan(thrust) # Pressure or thrust is NaN
        if tf1 or tf2 or tf3:
            break
    return df

"----- Inital Conditions -----"
"Injector"
Ainj = 3.4e-6                           # Area of injection holes | [m^2] 34mm from chinese
Cd = 0.5                                # Coefficient of Discharge, Injector | dimensionless

"Fuel Properties"
rho_fuel = 1166.15439765                # Density of fuel | [kg/m^3]
M = 28.811                              # Molecular Mass | [g / mol]
T = 5000                                # Chamber Temperature | [K]
gamma = 1.5                             # specific heat ratios
n = 0.8                                 # Pressure exponent | scalar
a = units(0.7, "ai_to_am", n=n)         # Burn rate coefficient | scalar  

"Combustion Chamber and Grain"
Fuel_OD = units(1.688, "in_to_m")       # Fuel Grain Outer Diameter | [in --> m] 
Fuel_ID = units(1.25, "in_to_m")        # Fuel Grain Inner Diameter | [in --> m]
L = units(12, "in_to_m")                # Fuel port length | [in --> m]
r_port0 = Fuel_OD - Fuel_ID             # Initial fuel port radius | [in --> m]
V_chmb_emty = units(120, "in3_to_m3")   # Empty Volume of the Combustion Chamber | [in^3 --> m^3]

"Oxidizer Tank"
m_loaded = 2.2                          # Initial mass of liquid oxidizer | kg
n_gas = 0.00028                         # Initial mass of Ullage gas | [kmol]
V_tank = units(3, "L_to_m3")            # Oxidizer Tank Volume | [L --> m^3]
Ti_tank = 298                           # Initial temperature of the oxidizer tank | K
m_T = 2.1;                              # Oxidizer tank mass [kg]

"Nozzle"
d_t = units(0.5, "in_to_m")             # Nozzle Throat diameter | [in --> m]
cd_throat = 0.2                         # Coefficient of Discharge of Nozzle Throat | dimensionless

"----- INPUTS -----"
data_columns = [ # Dataframe columns
    "time", 
    "dm_ox_dt", 
    "dr_dt", 
    "thrust",
    "n_lo",
    "n_go",
    "OF",
    "r",
    "V_chmb",
    "p_chmb",
    "p_tank"
]
df = pd.DataFrame(columns=data_columns) # Create dataframe

initial_inputs = { # Dictionary of initial inputs
    # Injector
    "Ainj": Ainj,
    "Cd": Cd,
    # Fuel Properties
    "rho_fuel": rho_fuel,
    "M": M,
    "T": T,
    "gamma": gamma,
    "n": n,
    "a": a,
    # Fuel Grain
    "L": L,
    "r_port0": r_port0,
    # Oxidizer Tank
    "m_loaded": m_loaded,
    "n_gas": n_gas,
    "V_tank": V_tank,
    "Ti_tank": Ti_tank,
    "m_T": m_T,
    # Nozzle
    "d_t": d_t,
    "cd_throat": cd_throat,
    # Volumes
    "V_chmb_emty": V_chmb_emty,
}

tf=20                           # final time [s]
tstep=0.005                     # time step [s]

df = mainLoop(tf, tstep, initial_inputs)

"----- Table -----"
print(tabulate(df.iloc[::10], headers='keys', tablefmt='fancy_grid')) # Print every 10th row
df.to_csv("src/Model/CSVFiles/current_data.csv", index=False) # Save dataframe to csv

"----- Plotting -----"
# plot_as_subplots(df) # Plot as subplots
plot_as_individual(df) # Plot as individual plots