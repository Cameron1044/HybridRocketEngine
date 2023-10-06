import numpy as np
import math as m
import matplotlib.pyplot as plt
from chemicalProperties import ChemicalProperties

"----- Inital Conditions -----"
"Oxidizer Tank - Liquid Nitroux Oxide"
R = 8314.3                              # universal gas constant [J/(kmol*K)]
MW2 = 44.013                            # molecular weight of N2O [kg / kmol]
rho_ox = 1226                           # Liquid Nitroux Oxide at 273K [kg/m³]
V_tank = 15.8213 / 61020                # Oxidizer Tank Volume | [in^3 --> m^3]

## DESIGN INPUTS:
V_ox_liquid = V_tank*(80/100)           # Initial volume of liquid oxidizer | [L or 1000 cm^3 --> m^3]
    # Note - assuming 80% of the volume is liquid nitrous oxide
m_loaded = 0.4/2.205 # V_ox_liquid * rho_ox   # Initial mass of liquid oxidizer | kg
V_gas  = V_tank*(20/100)                # Initial volume of Ullage Gas | [m^3]
n_gas = 0.000125                        # Initial mass of Ullage gas | [kmol]
P_ox = 2300 * 6895                      # Initial pressume of the oxidizer tank | [Psi --> Pa]
Ti_tank = 286.5                         # Initial temperature of the oxidizer tank | K
m_T = 6.4882;                           # Oxidizer tank mass [kg]

chem = ChemicalProperties(gas="He", mass_loaded=m_loaded, tank_volume=V_tank, initial_temp=Ti_tank)

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
a_m = a_i * conversion_factor_a           # Burn rate coefficient | [m/s]
## NOZZLE
d_t = 0.23 / 39.37                      # Nozzle Throat diameter | [in --> m]
cd_throat = 0.2                         # Coefficient of Discharge of Nozzle Throat | dimensionless


"----- Initial Calculations/Inputs -----"
Ainj = n_holes * np.pi * (phi/2)**2                             # Area of injection holes | [m^2]

# INITAL OXIDIZER TANK
n_go, n_lo = chem.initialMoles()
To = Ti_tank

# INITIAL COMBUSTION CHAMBER
A_t = np.pi*(d_t/2)**2                                          # Cross-sectional area of nozzle throat | [m^2]
A_port = 2 * np.pi * r_port0 * L                                # Inital fuel grain burn surface area | [m^2]
r = r_port0                                                     # Initial fuel grain radius | [m]

"----- CALCULATIONS -----"
# Time Loop
tf=10                           # final time [s]
tstep=0.031                     # time step [s]
i_f=tf/tstep

# Initial derivatives
dr_dt = 0
dm_ox_dt = 0
dV_dt = 0

# Preallocating arrays
t_arr = []
p_tank_arr = []
p_chmb_arr = []
thrust_arr = []
n_lo_arr = []
n_go_arr = []

# Initiating For Loop
for i in range(0,int(i_f)):

    "----- BLOWDOWN -----"
    Vhat_l, CVhat_He, CVhat_g, CVhat_l, delta_Hv, P_sat, dP_sat, Cp_T = chem.chemicalStates(To)

    ## Simplified expression definitions for solution
    P = (n_gas + n_go)*R*To / (V_tank - n_lo*Vhat_l)            # Calculating Oxidizer Tank Pressure
    a = m_T*Cp_T + n_gas*CVhat_He + n_go*CVhat_g + n_lo*CVhat_l
    b = P*Vhat_l
    e = -delta_Hv + R*To
    f = -Cd*Ainj*np.sqrt(2/MW2)*np.sqrt((P-P_chmb)/Vhat_l)
    j = -Vhat_l*P_sat
    k = (V_tank - n_lo*Vhat_l)*dP_sat
    m = R*To
    q = R*n_go

    Z=(-f*(-j*a + (q-k)*b)) / (a*(m+j) + (q-k)*(e-b))
    W=(-Z*(m*a + (q-k)*e)) / (-j*a + (q-k)*b)
    
    # Derivative Functions
    dT = (b*W+e*Z)/a
    dn_g = Z
    dn_l = W
    # Forward Difference Method
    To = To + dT*tstep
    n_go = n_go + dn_g*tstep
    n_lo = n_lo + dn_l*tstep
    # Saving Blowdown Model Results 
    n_lo_arr.append(n_lo)           # Moles of Liquid Nitrous Oxide
    n_go_arr.append(n_go)           # Moles of Gas Nitrous Oxide 
    p_tank_arr.append(P  / 6895)    # Pressure of the Oxidizer Tank [Pa --> Psi]

    "----- COMBUSTION CHAMBER -----"
    # Calculating Mass flow of Liquid Oxidizer
    dm_ox_dt = -1*dn_l*MW2                                              # [kg/s]
    m_ox = dm_ox_dt*tstep                                               # [kg]
    # Calculating change in fuel regression radius
    dr_dt = a_m * (dm_ox_dt / A_port)**n                                  # [m/s]
    r = r + dr_dt*tstep                                                 # [m]
    # Calculating Fuel Burn Surface Area
    A_port = 2 * np.pi * r * L                                          # [m^2]
    # Calculating mass of fuel flow rate
    m_fuel = dr_dt*tstep*A_port*rho_fuel                                # [kg]
    # Calculating Chamber Pressure
    C_star_nom = np.sqrt((R * T) / (gamma * M))                         # [m/s]
    C_star_denom = (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))
        # Note: This relation shows up in mdot exit and C_star | Dimensionless
    C_star = C_star_nom / C_star_denom  
    # Calculating Chamber Pressure
    P_chmb = ((dm_ox_dt + m_fuel/tstep) * C_star)/(A_t)                 # [Pa]
    # Calculating Chamber mass flow rate
    m_chmb = m_fuel + m_ox                                              # [kg/s]
    # Calculating new volume in the Combusation Chamber
    V_chmb = np.pi*((0.01587)**2)*(0.0127) + np.pi*((0.0158)**2)*(0.0254) + np.pi*((r)**2)*L #[in^3 --> m^3]
    # Calculating Density of the combustion chamber
    rho_chmb = m_chmb / V_chmb                                          # [kg/m^3]
    # Saving Combustion Chamber Results
    p_chmb_arr.append(P_chmb  / 6895) # Pressure of the Combustion Chamber [Pa --> Psi]

    "----- NOZZLE -----"
    # Calculating the choked mass flow of the combustion chamber
    CompressibleFactor = (2 / (gamma + 1))**((gamma + 1) / (gamma - 1))
    mdot_exit = cd_throat * A_t * np.sqrt(gamma * rho_chmb * P_chmb * CompressibleFactor) # [kg/s]
    # Calculating exit velocity
    v_exit = np.sqrt(((2*gamma)/(gamma-1)) * ((R*T)/M))
        # Note: Assumption for perfectly expanded nozzle such that P Chamber = P Exit
    # Calculating thrust
    thrust = mdot_exit * v_exit                                         # [N]
    # Saving Nozzle Results
    thrust_arr.append(thrust / 4.448)   # Thrust [N --> lbf]
    # Saving Time
    t_arr.append(i*tstep)
    # If chamber exceeds oxidizer tank pressure stop loop.
    if P < P_chmb:
        break
    # If mass of Liquid Nitroux Oxide is 0, stop loop
    if n_lo <= 0:
        break
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

# Oxidizer Tank Equillibrium of Nitroux Oxide Mass
plt.figure()
plt.grid(True)
plt.plot(t_arr, n_go_arr, 'b', linewidth=2)
plt.plot(t_arr, n_lo_arr, 'g', linewidth=2)
plt.title('Mass of N20 vs. Time')
plt.xlabel('Time [s]')
plt.ylabel('Mass of N2O [kg]')
plt.legend(['Mass of N2O gas', 'Mass of N2O liquid'])
plt.show()