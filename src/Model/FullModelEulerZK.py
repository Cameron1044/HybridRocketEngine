import numpy as np
import math as m
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
m_loaded = 0.4/2.205 # V_ox_liquid * rho_ox   # Initial mass of liquid oxidizer | kg
V_gas  = V_tank*(20/100)                # Initial volume of Ullage Gas | [m^3]
n_gas = 0.000125                        # Initial mass of Ullage gas | [kmol]
P_ox = 2300 * 6895                      # Initial pressume of the oxidizer tank | [Psi --> Pa]
Ti_tank = 286.5                         # Initial temperature of the oxidizer tank | K
m_T = 6.4882;                           # Oxidizer tank mass [kg]

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

"----- Chemical Properties -----"
# N20 Vapor
G1 = 96.512                             # vapor pressure of N2O [Pa] coefficients
G2 = -4045                              # valid for Temp range [182.3 K - 309.57 K]
G3 = -12.277
G4 = 2.886e-5
G5 = 2

D1 = 0.2934e5                           # heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients
D2 = 0.3236e5                           # valid for Temp range [100 K - 1500 K]
D3 = 1.1238e3
D4 = 0.2177e5
D5 = 479.4

Tc = 309.57                             # critical temperature of N2O [K]
J1 = 2.3215e7                           # heat of vaporization of N2O [J/kmol] coefficients
J2 = 0.384                              # valid for Temp range [182.3 K - 309.57 K]
J3 = 0
J4 = 0
# Helium
C1 = 0.2079e5                           # heat capacity of He at constant pressure [J/(kmol*K)] coefficients
C2 = 0                                  # valid for Temp range [100 K - 1500 K]
C3 = 0                                  # Same heat capacity as Argon Page 182 or PDFPage 207
C4 = 0
C5 = 0
# N20 Liquid
Q1 = 2.781                              # molar specific volume of liquid N2O [m**3/kmol] coefficients
Q2 = 0.27244
Q3 = 309.57
Q4 = 0.2882
E1 = 6.7556e4                           # heat capacity of N2O liquid at constant pressure [J/(kmol*K)] coefficients
E2 = 5.4373e1                           # valid for Temp range [182.3 K - 200 K]
E3 = 0
E4 = 0
E5 = 0

"----- Initial Calculations/Inputs -----"
Ainj = n_holes * np.pi * (phi/2)**2                             # Area of injection holes | [m^2]

# INITAL OXIDIZER TANK
n_to = m_loaded / MW2                                           # initial total N2O in tank [kmol]
Vhat_li = Q2**(1 + (1 - Ti_tank / Q3)**Q4) / Q1                 # molar volume of liquid N2O [m**3/kmol]
To = Ti_tank                                                    # initial temperature [K]
P_sato = np.exp(G1 + G2 / To + G3 * np.log(To) + G4 * To**G5)   # initial vapor pressure of N20 [Pa]
n_go = P_sato * (V_tank - Vhat_li * n_to) / (-P_sato * Vhat_li + R * To)     # initial N2O gas [kmol]
n_lo = (n_to * R * To - P_sato * V_tank) / (-P_sato * Vhat_li + R * To)      # initial N2O liquid [kmol]


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
    # Given functions of temperature:
    Vhat_l = Q2**(1+(1-To/Q3)**Q4)/Q1                                                                   # molar specific volume of liquid N2O [m^3/kmol]

    ## Calculating Specific Heats at Constant Volume
    CVhat_He = C1 + C2*To + C3*To**2 + C4*To**3 + C5*To**4 - R                                          # specific heat of He [J/(kmol*K)]
    CVhat_g = D1 + D2*((D3/To)/np.sinh(D3/To))**2 + D4*((D5/To)/np.cosh(D5/To))**2 - R                  # specific heat of N2O gas [J/(kmol*K)]
    CVhat_l = E1 + E2*To + E3*To**2 + E4*To**3 + E5*To**4                                               # specific heat of N2O liquid [J/(kmol*K)]
        # Note: At constant volume, approx. same as at constant pressure 

    Tr = To/Tc # reduced temperature
    delta_Hv = J1*(1 - Tr) ** (J2 + J3*Tr + J4*Tr**2)                                                   # heat of vaporization of N2O [J/kmol]
    P_sat = np.exp(G1 + G2/To + G3*np.log(To) + G4*To**G5)                                              # vapor pressure of N20 [Pa]
    dP_sat = (-G2/(To**2) + G3/To + G4*G5*To**(G5-1)) * np.exp(G1 + G2/To + G3*np.log(To) + G4*To**G5)  # derivative of vapor pressure with respect to temperature
    Cp_T = (4.8 + 0.00322*To)*155.239                                                                   # specific heat of oxidizer tank, Aluminum [J/(kg*K)]

    ## Simplified expression definitions for solution
    P = (n_gas + n_go)*R*To / (V_tank - n_lo*Vhat_l)                                                    # Calculating Oxidizer Tank Pressure
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
    print(P / 6895)
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
    print(dm_ox_dt)
    dr_dt = a * (dm_ox_dt / A_port)**n                                  # [m/s]
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
    print((P_chmb  / 6895))
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