import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

"----- Inital Conditions -----"
"Oxidizer Tank - Liquid Nitroux Oxide"
R = 8314.3                              # universal gas constant [J/(kmol*K)]
MW2 = 44.013                            # molecular weight of N2O [kg / kmol]
rho_ox = 1226                           # Liquid Nitroux Oxide at 273K [kg/m³]
V_tank = 0.003                          # Oxidizer Tank Volume | [in^3 --> m^3]

## DESIGN INPUTS:
V_ox_liquid = V_tank*(80/100)           # Initial volume of liquid oxidizer | [L or 1000 cm^3 --> m^3]
    # Note - assuming 80% of the volume is liquid nitrous oxide
P_ox = 2300 * 6895                      # Initial pressume of the oxidizer tank | [Psi --> Pa]

"Injector"
Cd = 0.5                                # Coefficient of Discharge, Injector | dimensionless
## Commented to do later
# n_holes = 4                             # Number of holes, injector | quantity
# phi = 0.0492 / 39.37                    # Diameter of injection holes | [in --> m]

"Combustion Chamber Fuel Grain"
P_chmb = 0                              # Initial Pressure of Combustion Chamber | [Psi --> Pa]
V_chmb_emty = 120 / 61020               # Empty Volume of the Combustion Chamber | [in^3 --> m^3]
## DESIGN INPUTS OF FUEL:
L = 12 / 39.37                          # Fuel port length | [in --> m]
Fuel_OD = 1.688 / 39.37                 # Fuel Grain Outer Diameter | [in --> m] 
Fuel_ID = 1.25 / 39.37                  # Fuel Grain Inner Diameter | [in --> m]   
r_port0 = Fuel_OD - Fuel_ID             # Initial fuel port radius | [in --> m]
## Burned Oxidizer + Fuel Properties
gamma = 1.5                             # specific heat ratios
T = 5000                                # Chamber Temperature | [K]
M = 28.811                              # Molecular Mass | [g / mol]
rho_fuel = 1166.15439765                # Density of fuel | [kg/m^3]
# Fuel Regression Coefficients:
a_i = 0.7                               # Burn rate coefficient | scalar     
n = 0.8                                 # Pressure exponent | scalar
conversion_factor_a = (0.0254**(1+2*(n)))*(0.453592**(-n)) # Conversion factor for a | [in --> m, lb --> kg]
a_m = a_i * conversion_factor_a           # Burn rate coefficient | [m/s]
## NOZZLE
d_t = 0.013                             # Nozzle Throat diameter | [in --> m]
cd_throat = 0.2                         # Coefficient of Discharge of Nozzle Throat | dimensionless

"----- Initial Calculations/Inputs -----"
Ainj = 3.4e-6                           # Area of injection holes | [m^2] 34mm
    # n_holes * np.pi * (phi/2)**2

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
tf=15                           # final time [s]
tstep=0.005                     # time step [s]
i_f=tf/tstep

# Initial derivatives
dr_dt = 0
dm_ox_dt = mdot_ox0
dV_dt = 0

# Preallocating arrays
t_arr = []              # time
p_tank_arr = []         # Oxidizer tank pressure
p_chmb_arr = []         # Combustion chamber pressure
thrust_arr = []         # Thrust
OF_arr = []             # Oxidizer to Fuel Ratio
r_arr = []              # Fuel Grain Radius
V_chmb_arr = []         # Combustion Chamber Volume
mo_arr = []             # Mass flow of Oxidizer
mf_arr = []             # Mass flow of Fuel

# Initiating For Loop
for i in range(0,int(i_f)):
    # Calculating Pressure of Oxidizer Tank
    P_ox = P1V1 / V_ox_gas

    # If chamber exceeds oxidizer tank pressure stop loop.
    if P_ox < P_chmb:
        break
    # Calculating Mass flow of Liquid Oxidizer
    dm_ox_dt = Cd * Ainj * np.sqrt(2 * rho_ox * (P_ox - P_chmb))        # [kg/s]
    m_ox = dm_ox_dt*tstep                                               # [kg]
    # Calculating change in Volume of Oxidizer Tank
    dV_dt = dm_ox_dt / rho_ox * tstep                                   # [m^3]
    V_ox_gas = V_ox_gas + dV_dt                                         # [m^3]
    # Calculating change in fuel regression radius
    dr_dt = a_m * (dm_ox_dt / A_port)**n                                # [m/s]
    r = r - dr_dt*tstep                                                 # [m]
    # Calculating Fuel Burn Surface Area
    A_port = 2 * np.pi * r * L                                          # [m^2]
    # Calculating mass of fuel flow rate
    dm_f_dt = dr_dt * A_port * rho_fuel                                 # [kg/s] 
    m_fuel = dm_f_dt * tstep                                            # [kg]
    # Calculating Oxidizer to Fuel Ratio
    OF = dm_ox_dt / dm_f_dt
    # Calculating Chamber Pressure
    C_star_nom = np.sqrt((R * T) / (gamma * M))                         # [m/s]
    C_star_denom = (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))
        # Note: This relation shows up in mdot exit and C_star | Dimensionless
    C_star = C_star_nom / C_star_denom  
    # Calculating Chamber Pressure
    p_chmb = ((dm_ox_dt + m_fuel/tstep) * C_star)/(A_t)                 # [Pa]
    # Saving Combustion Chamber Results
    p_tank_arr.append(P_ox  / 6895)    # Pressure of the Oxidizer Tank [Pa --> Psi]
    r_arr.append(r * 39.37)            # Fuel Grain Radius [m --> in]
    OF_arr.append(OF)                  # Fuel to oxidizer ratio
    mo_arr.append(dm_ox_dt * 2.205)    # Oxidizer Mass flow rate [kg/s --> lb/s]
    mf_arr.append(dm_f_dt * 2.205)     # Fuel Mass flow rate [kg/s --> lb/s]
    p_chmb_arr.append(p_chmb  / 6895)  # Pressure of the Combustion Chamber [Pa --> Psi]

    # Calculating Chamber mass flow rate
    m_chmb = m_fuel + m_ox                                              # [kg/s]
    # Calculating new volume in the Combusation Chamber
    V_chmb = V_chmb_emty - (np.pi * r**2 * L)                           # [m^3]
    # Calculating Density of the combustion chamber
    rho_chmb = m_chmb / V_chmb                                          # [kg/m^3]
    # Calculating the choked mass flow of the combustion chamber
    CompressibleFactor = (2 / (gamma + 1))**((gamma + 1) / (gamma - 1))
    mdot_exit = cd_throat * A_t * np.sqrt(gamma * rho_chmb * p_chmb * CompressibleFactor) # [kg/s]
    # Calculating exit velocity
    v_exit = np.sqrt(((2*gamma)/(gamma-1)) * ((R*T)/M))
        # Note: Assumption for perfectly expanded nozzle such that P Chamber = P Exit
    # Calculating thrust
    thrust = mdot_exit * v_exit                                         # [N]
    # Saving Nozzle Results
    V_chmb_arr.append(V_chmb * 61020)  # Volume of the Combustion Chamber [m^3 --> in^3]
    thrust_arr.append(thrust / 4.448)   # Thrust [N --> lbf]
    # Saving Time
    t_arr.append(i*tstep)
# End of For Loop

"----- Plotting -----"
# Thrust Profile
plt.figure()
plt.grid(True)
plt.plot(t_arr, thrust_arr, linewidth=2)
plt.title('Thrust vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Thrust (lbf)')
plt.show()

# Oxidizer Tank and Combustion Chamber Pressure
plt.figure()
plt.grid(True)
plt.plot(t_arr, p_chmb_arr, linewidth=2, label="Combustion Chamber")
plt.plot(t_arr, p_tank_arr, linewidth=2, label="Oxidizer Tank") 
plt.title('Oxidizer Tank and Combustion Chamber Pressure vs. Time')
plt.legend(loc = 'best')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.show()

# OF Profile
plt.figure()
plt.grid(True)
plt.plot(t_arr, OF_arr, linewidth=2)
plt.title('Oxidizer to Fuel Ratio vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Oxidizer to Fuel Ratio (O/F)')
plt.show()

# Fuel Grain Radius Profile
plt.figure()
plt.grid(True)
plt.plot(t_arr, r_arr, linewidth=2)
plt.title('Fuel Grain Burn over Time')
plt.xlabel('Time (s)')
plt.ylabel('Fuel Grain Burn Rate (r) [in]')
plt.show()

# Combustion Chamber Profile
plt.figure()
plt.grid(True)
plt.plot(t_arr, V_chmb_arr, linewidth=2)
plt.title('Combustion Chamber Volume vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Combustion Chamber Volume [in^3]')
plt.show()

# Mass Flow Rate of Oxidizer and Fuel
plt.figure()
plt.grid(True)
plt.plot(t_arr, mo_arr, linewidth=2, label='Oxidizer')
plt.plot(t_arr, mf_arr, linewidth=2, label='Fuel Grain')
plt.title('Mass Flow Rates over Time')
plt.legend(loc='best')
plt.xlabel('Time (s)')
plt.ylabel('Mass Flow Rates [lb/s]')
plt.show()