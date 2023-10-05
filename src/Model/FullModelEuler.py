import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

"----- Inital Conditions -----"
"Oxidizer Tank - Liquid Nitroux Oxide"
R = 8314.3                              # universal gas constant [J/(kmol*K)]
MW2 = 44.013                            # molecular weight of N2O [kg / kmol]
rho_ox = 1226                           # Liquid Nitroux Oxide at 273K [kg/m³]
V_tank = 15.8213 / 61020                # Oxidizer Tank Volume | [in^3 --> m^3]

## DESIGN INPUTS:
V_ox_liquid = V_tank*(80/100)           # Initial volume of liquid oxidizer | [L or 1000 cm^3 --> m^3]
    # Note - assuming 80% of the volume is liquid nitrous oxide
P_ox = 2300 * 6895                      # Initial pressume of the oxidizer tank | [Psi --> Pa]

"Injector"
Cd = 0.175                              # Coefficient of Discharge, Injector | dimensionless
n_holes = 4                             # Number of holes, injector | quantity
phi = 0.0492 / 39.37                    # Diameter of injection holes | [in --> m]

"Combustion Chamber Fuel Grain"
P_chmb = 0                              # Initial Pressure of Combustion Chamber | [Psi --> Pa]
## DESIGN INPUTS OF FUEL:
L = 4 / 39.37                           # Fuel port length | [in --> m]
r_port0 = 0.25 / 39.37                  # Initial fuel port radius | [in --> m]
## Burned Oxidizer + Fuel Properties
gamma = 1.249                           # specific heat ratios
T = 3059                                # Chamber Temperature | [K]
M = 28.811                              # Molecular Mass | [g / mol]
rho_fuel = 1166.15439765                # Density of fuel | [kg/m^3]
# Fuel Regression Coefficients:
a_i = 0.7                               # Burn rate coefficient | scalar     
n = 0.8                                 # Pressure exponent | scalar
conversion_factor_a = (0.0254**(1+2*(n)))*(0.453592**(-n)) # Conversion factor for a | [in --> m, lb --> kg]
a = a_i * conversion_factor_a           # Burn rate coefficient | [m/s]
## NOZZLE
d_t = 0.23 / 39.37                      # Nozzle Throat diameter | [in --> m]
cd_throat = 0.2                         # Coefficient of Discharge of Nozzle Throat | dimensionless

"----- Initial Calculations/Inputs -----"
Ainj = n_holes * np.pi * (phi/2)**2                             # Area of injection holes | [m^2]

# INITAL OXIDIZER TANK
mdot_ox0 = Cd*Ainj*np.sqrt(2*rho_ox*(P_ox - P_chmb))             # Initial mass flow rate of oxidizer | [kg/s]
V_ox_gas0 = V_tank - V_ox_liquid                                 # Initial volume of gas pressurant | [m^3]
P1V1 = P_ox * V_ox_gas0                                          # Initial pressure * volume in oxidizer tank

# INITIAL COMBUSTION CHAMBER
A_t = np.pi*(d_t/2)**2                                          # Cross-sectional area of nozzle throat | [m^2]
A_port = 2 * np.pi * r_port0 * L                                # Inital fuel grain burn surface area | [m^2]
r = r_port0                                                     # Initial fuel grain radius | [m]
V_ox_gas = V_ox_gas0                                            # Initial Volume of Oxidizer Tank Gas

"----- CALCULATIONS -----"
# Time Loop
tf=10                           # final time [s]
tstep=0.031                     # time step [s]
i_f=tf/tstep

# Initial derivatives
dr_dt = 0
dm_ox_dt = mdot_ox0
dV_dt = 0

# Preallocating arrays
t_arr = []
p_chmb_arr = []
thrust_arr = []

# Initiating For Loop
for i in range(0,int(i_f)):
    # Calculating Pressure of Oxidizer Tank
    P_ox = P1V1 / V_ox_gas

    # If chamber exceeds oxidizer tank pressure stop loop.
    if P_ox < P_chmb:
        break

    # If gas volume takes up the volume of the oxidizer tank
    # if V_tank < V_ox_gas:
    #     break 

    # Calculating Mass flow of Liquid Oxidizer
    dm_ox_dt = Cd * Ainj * np.sqrt(2 * rho_ox * (P_ox - P_chmb))        # [kg/s]
    m_ox = dm_ox_dt*tstep                                               # [kg]
    # Calculating change in Volume of Oxidizer Tank
    dV_dt = dm_ox_dt / rho_ox * tstep                                   # [m^3]
    V_ox_gas = V_ox_gas + dV_dt                                         # [m^3]
    # Calculating change in fuel regression radius
    dr_dt = a * (dm_ox_dt / A_port)**n                                  # [m/s]
    r = r + dr_dt*tstep                                                 # [m]
    # Calculating Fuel Burn Surface Area
    A_port = 2 * np.pi * r * L                                          # [m^2]
    # Calculating mass of fuel flow rate
    m_fuel = dr_dt*tstep*A_port*rho_fuel                                # [kg]
    # Calculating Chamber Pressure
    C_star_nom = np.sqrt((R * T) / (gamma * M))                   # [m/s]
    CompressibleFactor = (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))
        # Note: This relation shows up in mdot exit and C_star | Dimensionless
    C_star = C_star_nom / CompressibleFactor  

    p_chmb = ((dm_ox_dt + m_fuel/tstep) * C_star)/(A_t)                 # [Pa]
    # p_chmb = ((dm_ox_dt + m_fuel/tstep))*(840.8786447/0.588693023)/(A_t)
    # print(p_chmb1, p_chmb, C_star_nom, CompressibleFactor, C_star)

    t_arr.append(i*tstep)
    p_chmb_arr.append(p_chmb  / 6895) # Converting Pa to Psi

    # Calculating Chamber mass flow rate
    m_chmb = m_fuel + m_ox                                              # [kg/s]
    # Calculating new volume in the Combusation Chamber
    V_chmb = np.pi*((0.01587)**2)*(0.0127) + np.pi*((0.0158)**2)*(0.0254) + np.pi*((r)**2)*L #[in^3 --> m^3]
    # print(V_chmb * 61020)
    # Calculating Density of the combustion chamber
    rho_chmb = m_chmb / V_chmb                                          # [kg/m^3]
    # Calculating the choked mass flow of the combustion chamber
    mdot_exit = cd_throat * A_t * np.sqrt(gamma * rho_chmb * p_chmb * CompressibleFactor) # [kg/s]
    # Calculating exit velocity
    v_exit = np.sqrt(((2*gamma)/(gamma-1)) * ((R*T)/M))
        # Note: Assumption for perfectly expanded nozzle such that P Chamber = P Exit
    # Calculating thrust
    thrust = mdot_exit * v_exit                                         # [N]
    thrust_arr.append(thrust / 4.448) # Converting N to lbf
    print(i*tstep, dm_ox_dt*2.20462, dr_dt*39.37, r*39.37, m_fuel*2.20462, p_chmb/6895)

    #print(dm_ox_dt, dV_dt, r, p_chmb, thrust)
# End of For Loop

"----- Plotting -----"
# Chamber Pressure
plt.figure()
plt.grid(True)
plt.plot(t_arr, p_chmb_arr) 
plt.xlabel('Time (s)')
plt.ylabel('P_chmb (psi)')
plt.show()

# Thrust Profile
plt.figure()
plt.grid(True)
plt.plot(t_arr, thrust_arr)
plt.xlabel('Time (s)')
plt.ylabel('Thrust (lbf)')
plt.show()