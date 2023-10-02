import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

"----- Inital Conditions -----"
"Oxidizer Tan - Liquid Nitroux Oxide"
rho_ox = 0.0440                     # Density of oxidizer | lb/in^3
p_ox = 2300                         # Pressure of oxidizer tank | psi
V_ox = 15.8213                      # Total volume of oxidizer tank | in^3
V_ox_liquid = V_ox*(80/100)         # Initial volume of liquid oxidizer | in^3
    # Note - assuming 80% of the volume is liquid nitrous oxide

"Injector"
Cd = 0.175                          # Coefficient of discharge, orifice | dimensionless
n_holes = 4                         # Number of holes, injector | quantity
phi = 0.0492                        # Diameter of injection holes | inches

"Combustion Chamer Fuel Grain"
p_chmb = 0                          # Pressure of Combustion Chamber, pre-fire | psi
L = 4                               # Fuel port length | inches
r_port0 = 0.25                      # Initial fuel port radius | inches
a = 0.7                             # Burn rate coefficient | scalar
n = 0.8                             # Pressure exponent | scalar
d_t = 0.23                          # Throat diameter | inches
rho_fuel = 0.0421                   # Density of fuel | lb/in^3

### Constants ###
gamma = 1.249                       # specific heat ratios
cd_throat = 0.2
Ru = 8.314                          # Universal Gas COnstanty
T = 3059
M = 28.811                          # Molecular Mass

"----- Initial Calculations/Inputs -----"
A = n_holes*np.pi*(phi/2)**2                                    # Area of injection hole | in^2

# INITAL OXIDIZER TANK
mdot_ox0 = Cd*A*np.sqrt(2*rho_ox*(p_ox - p_chmb)*32.2*12)       # Initial mass flow rate of oxidizer | lb/s
V_ox_gas0 = V_ox - V_ox_liquid                                  # Initial volume of gas pressurant | in^3
P1V1 = p_ox * V_ox_gas0                                         # Initial pressure * volume in oxidizer tank

# INITIAL COMBUSTION CHAMBER
A_t = np.pi*(d_t/2)**2                                          # Cross-sectional area of nozzle throat | in^2
A_port = 2 * np.pi * r_port0 * L                                # Inital fuel grain burn surface area | in^2

r = r_port0
V_ox_gas = V_ox_gas0

"----- CALCULATIONS -----"
# Time Step
tf=100
tstep=0.031
i_f=tf/tstep

# Initial derivatives
dr_dt = 0
dm_ox_dt = mdot_ox0
dV_dt = 0

# Preallocating arrays
t_arr = []
p_chmb_arr = []
thrust_choked_arr = []

# Initiating For Loop
for i in range(0,int(i_f)):
    # Calculating Pressure of Oxidizer Tank
    p_ox = P1V1 / V_ox_gas

    # If chamber exceeds oxidizer tank pressure stop loop.
    if p_ox < p_chmb:
        break

    # Calculating Mass flow of Liquid Oxidizer
    dm_ox_dt = Cd * A * np.sqrt(2 * rho_ox * (p_ox - p_chmb) * 32.2 * 12)
    m_ox = dm_ox_dt*tstep
    # Calculating change in Volume
    dV_dt = dm_ox_dt / rho_ox * tstep
    V_ox_gas = V_ox_gas + dV_dt
    # Calculating change in fuel regression radius
    dr_dt = a * (dm_ox_dt / A_port)**n 
    r = r + dr_dt*tstep
    # Calculating Fuel Burn Surface Area
    A_port = 2 * np.pi * r * L
    # Calculating mass of fuel flow rate
    m_fuel = dr_dt*tstep*A_port*rho_fuel
    # Calculating Chamber Pressure
    p_chmb = (3.28084*(dm_ox_dt + m_fuel/tstep))*(840.8786447/0.588693023)/(32.2*A_t)

    t_arr.append(i*tstep)
    p_chmb_arr.append(p_chmb)

    # Calculating Chamber mass flow rate
    m_chmb = m_fuel + m_ox
    # Calculating open volume in the Combusation Chamber
    V_chmb = np.pi*(1.25/2)**2*0.5 + np.pi*(1.25/2)**2*1 + np.pi*(r)**2*L
    # Calculating Density of the combustion chamber
    rho_chmb = m_chmb / V_chmb
    # Calculating the chocked mass flow of the combustion chamber
    mdot_exit_choked = cd_throat * A_t * np.sqrt(32.2*12*gamma*rho_chmb*p_chmb*(2/(gamma+1))**((gamma+1)/(gamma-1)))
    # Calculating exit velocity
    v_exit = 3.28084*np.sqrt(((1000*2*gamma)/(gamma-1))*((Ru*T)/M))
    # Calculating thrust
    thrust_choked = mdot_exit_choked * v_exit / 32.2
    thrust_choked_arr.append(thrust_choked)

    print(dm_ox_dt, dV_dt, r, p_chmb, thrust_choked)
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
plt.plot(t_arr, thrust_choked_arr)
plt.xlabel('Time (s)')
plt.ylabel('Thrust (psi)')
plt.show()