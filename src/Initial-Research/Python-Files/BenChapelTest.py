import numpy as np
import math as m
import matplotlib.pyplot as plt
from chemicalProperties import ChemicalProperties

"----- Inital Conditions -----"
"Oxidizer Tank - Liquid Nitroux Oxide"
R = 8314.3                              # universal gas constant [J/(kmol*K)]
MW2 = 44.013                            # molecular weight of N2O [kg / kmol]
rho_ox = 1226                           # Liquid Nitroux Oxide at 273K [kg/m³]
V_tank = 0.003                          # Oxidizer Tank Volume | [in^3 --> m^3] 3L

# Meh
P_amb = 14 * 6895                       # Ambient Pressure [psi --> Pa]

## DESIGN INPUTS:
ullage_frac = 0.8
V_ox_liquid = V_tank*ullage_frac        # Initial volume of liquid oxidizer | [L or 1000 cm^3 --> m^3]
    # Note - assuming 80% of the volume is liquid nitrous oxide
m_loaded = V_ox_liquid * rho_ox         # Initial mass of liquid oxidizer | kg
V_gas  = V_tank*(1 - ullage_frac)       # Initial volume of Ullage Gas | [m^3]
n_gas = 0.0000125                       # Initial amount of Ullage gas | [kmol]
P_ox = 3000 * 6895                      # Initial pressume of the oxidizer tank | [Psi --> Pa]
Ti_tank = 298                           # Initial temperature of the oxidizer tank | K
m_T = 2.1;                              # Oxidizer tank mass [kg]

chem = ChemicalProperties(gas="Ar", mass_loaded=m_loaded, tank_volume=V_tank, initial_temp=Ti_tank)

"Injector"
Cd = 0.5                                # Coefficient of Discharge, Injector | dimensionless
## Commented to do later
# n_holes = 4                             # Number of holes, injector | quantity
# phi = 0.0492 / 39.37                    # Diameter of injection holes | [in --> m]

"Combustion Chamber Fuel Grain"
## DESIGN INPUTS OF FUEL:
L = 12 / 39.37                          # Fuel port length | [in --> m]
Fuel_OD = 3.376 / 39.37                 # Fuel Grain Outer Diameter | [in --> m] 
Fuel_ID = 1.25 / 39.37                  # Fuel Grain Inner Diameter | [in --> m]   
r_port0 = Fuel_ID / 2                   # Initial fuel port radius | [in --> m]
P_chmb = 0                              # Initial Pressure of Combustion Chamber | [Psi --> Pa]
V_chmb_emty = (Fuel_OD/2)**2 * np.pi * L              # Empty Volume of the Combustion Chamber | [in^3 --> m^3]
print(V_chmb_emty)
## Burned Oxidizer + Fuel Properties
gamma = 1.5                             # specific heat ratios
T = 3000                                # Chamber Temperature | [K]
M = 28.811                              # Molecular Mass | [g / mol]
rho_fuel = 1166.15439765                # Density of fuel | [kg/m^3]
# Fuel Regression Coefficients:
a_i = 0.7                               # Burn rate coefficient | scalar     
n = 0.8                                 # Pressure exponent | scalar
conversion_factor_a = (0.0254**(1+2*(n)))*(0.453592**(-n)) # Conversion factor for a | [in --> m, lb --> kg]
a_m = a_i * conversion_factor_a           # Burn rate coefficient | [m/s]
## NOZZLE
d_t = 0.013                             # Nozzle Throat diameter | [in --> m]
cd_throat = 0.85                         # Coefficient of Discharge of Nozzle Throat | dimensionless


"----- Initial Calculations/Inputs -----"
Ainj = 3.4e-6                           # Area of injection holes | [m^2] 34mm from chinese
    # n_holes * np.pi * (phi/2)**2
# INITAL OXIDIZER TANK
n_go, n_lo = chem.initialMoles()
To = Ti_tank

# INITIAL COMBUSTION CHAMBER
A_t = np.pi*(d_t/2)**2                                          # Cross-sectional area of nozzle throat | [m^2]
A_port = 2 * np.pi * r_port0 * L                                # Inital fuel grain burn surface area | [m^2]
r_port = r_port0                                                     # Initial fuel grain radius | [m]

"----- CALCULATIONS -----"
# Time Loop
tf=20                           # final time [s]
tstep=0.005                     # time step [s]
i_f=tf/tstep

# Initial derivatives
dr_dt = 0
dm_ox_dt = 0
dV_dt = 0

# Preallocating arrays
t_arr = [] # time
p_tank_arr = [] # Oxidizer tank pressure
p_chmb_arr = [] # Combustion chamber pressure
thrust_arr = [] # Thrust
n_lo_arr = [] # Moles of Liquid Nitrous Oxide
n_go_arr = [] # Moles of Gaseous Nitrous Oxide
OF_arr = [] # Oxidizer to Fuel Ratio
r_arr = [] # Fuel Grain Radius
V_chmb_arr = [] # Combustion Chamber Volume
mo_arr = [] # Mass

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
    """# Calculating Chamber Pressure
    C_star_nom = np.sqrt((R * T) / (gamma * M))                         # [m/s]
    C_star_denom = (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))
        # Note: This relation shows up in mdot exit and C_star | Dimensionless
    C_star = C_star_nom / C_star_denom  """
    # Calculating Chamber Pressure
    print(A_port, A_t, a_m, rho_fuel, R, T, n)
    P_chmb = (A_port/A_t * (a_m * rho_fuel)/np.sqrt(gamma/(R/M*T) * (2/(gamma+1))**((gamma+1)/(gamma-1))))**(1/(1-n))                    # [Pa]
    print(P_chmb,P)
    # Calculating Mass flow of Liquid Oxidizer
    dm_ox_dt = -1*dn_l*MW2                                              # [kg/s]
    print(dm_ox_dt)
    # Calculating change in fuel regression radius
    r = a_m * (P_chmb)**n                                # [m/s]
    r_port = r_port + r*tstep                                                 # [m]
    # Calculating Fuel Burn Surface Area
        # Note: Cylindrical Surface area --> 2 pi r h
    A_port = 2 * np.pi * r_port * L                                          # [m^2]
    # Calculating mass of fuel flow rate
    dm_f_dt = r * A_port * rho_fuel                                 # [kg/s]
    # Calculating Oxidizer to Fuel Ratio
    OF = dm_ox_dt / dm_f_dt
    # Calculating total mass flow rate
    dm_dt = dm_ox_dt + dm_f_dt                                               # [kg/s]
    # Calculating mass of combusting propellant that pases through chamber during timestep
    m_chmb = dm_dt * tstep                                         # [kg/s]
    # Calculating new volume in the Combusation Chamber
    V_chmb = V_chmb_emty - (np.pi * r_port**2 * L)                           # [m^3]
    # Calculating Density of the combustion chamber
    rho_chmb = m_chmb / V_chmb                                          # [kg/m^3]
    # Saving Combustion Chamber Results
    r_arr.append(r * 39.37)           # Fuel Grain Radius [m --> in]
    V_chmb_arr.append(V_chmb * 61020) # Volume of the Combustion Chamber [m^3 --> in^3]
    OF_arr.append(OF)                 # Fuel to oxidizer ratio
    p_chmb_arr.append(P_chmb  / 6895) # Pressure of the Combustion Chamber [Pa --> Psi]

    "----- NOZZLE -----"
    # Calculating the choked mass flow of the combustion chamber
     # CompressibleFactor = (2 / (gamma + 1))**((gamma + 1) / (gamma - 1))
     # mdot_exit = cd_throat * A_t * np.sqrt(gamma * rho_chmb * P_chmb * CompressibleFactor) # [kg/s]
    mdot_exit = dm_dt
    # Calculating exit velocity
    v_exit = np.sqrt(2 * T * (R/M) * (gamma/(gamma-1)) * (1 - (P_amb/P_chmb)**((gamma-1)/gamma)))  # [m/s] assuming optimal expansion
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

print(p_chmb_arr)

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
plt.plot(t_arr, p_chmb_arr, linewidth=2)
plt.plot(t_arr, p_tank_arr, linewidth=2) 
plt.title('Oxidizer Tank and Combustion Chamber Pressure vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.show()

# Oxidizer Tank Equillibrium of Nitroux Oxide Mass
plt.figure()
plt.grid(True)
plt.plot(t_arr, n_go_arr, 'b', linewidth=2)
plt.plot(t_arr, n_lo_arr, 'g', linewidth=2)
plt.title('Mass of N20 vs. Time')
plt.xlabel('Time [s]')
plt.ylabel('Mass of N2O [Moles]')
plt.legend(['Mass of N2O gas', 'Mass of N2O liquid'])
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