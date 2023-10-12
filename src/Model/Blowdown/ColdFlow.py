## Converted IdealLiquidModel.m to Python
## One diemnsional, unsteady compressible flow
import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

" DIMENSIONS OF ENGINE "
## Tank
m_loaded = 19.32933                     # N2O mass initially loaded into oxidizer tank [kg]: Test 1
V = 0.0354                              # total oxidizer tank volume [m**3]
## Combustion Chamber
A = 0.1                                 # Area of Injector Holes
V_chmb = 0.0400                         # total combustion chamber volume [m**3]
## Nozzle
d_t = 0.23                              # Throat diameter []

" OXIDIZER TANK CHEMICAL PROPERTIES "
R = 8314.3                              # universal gas constant [J/(kmol*K)]
MW2 = 44.013                            # molecular weight of N2O [g / mol]
rho_ox = 1.9781                         # Liquid Nitroux Oxide at 273K [kg/m³]
# Perry's Chemical Engineers' Handbook Property Equations
G1 = 96.512                             # vapor pressure of N2O [Pa] coefficients
G2 = -4045                              # valid for Temp range [182.3 K - 309.57 K]
G3 = -12.277
G4 = 2.886e-5
G5 = 2
# molar specific volume of liquid N2O [m**3/kmol] coefficients
Q1 = 2.781                              
Q2 = 0.27244
Q3 = 309.57
Q4 = 0.2882

" Fuel Grain "
L = 25                                  # Fuel Port Length []
r_port0 = 0.25                          # Initial fuel port radius []
a = 0.7                                 # Burn rate coefficient [scalar]
n = 0.8                                 # Pressure exponent [scalar]
rho_fuel = 0.0421                       # Density of Fuel Grain []

"-------- START VALUES ---------"
Ti = 286.5                              # initial temperature [K]: Test 1

### Calculating Initial Conditions of Oxidizer Tank ###
n_to = m_loaded / MW2                                                   # initial total N2O in tank [kmol]
Vhat_li = Q2**(1 + (1 - Ti / Q3)**Q4) / Q1                              # molar volume of liquid N2O [m**3/kmol]
To = Ti                                                                 # initial temperature [K]
P_sato = np.exp(G1 + G2 / To + G3 * np.log(To) + G4 * To**G5)           # initial vapor pressure of N20 [Pa]
n_go = P_sato * (V - Vhat_li * n_to) / (-P_sato * Vhat_li + R * To)     # initial N2O gas [kmol]
n_lo = (n_to * R * To - P_sato * V) / (-P_sato * Vhat_li + R * To)      # initial N2O liquid [kmol]

### Calculating Initial Conditions of Combustion Chamber ###
A_port = 2 * np.pi * r_port0 * L                                        # Initial Fuel Grain exposed area
r = r_port0                                                             # Initial radius of Fuel Grain

# Initial Conditions passed into ode45 state vector
y0 = [Ti, n_go, n_lo, r]

"-------- FUNCTIONS ---------"
### Defining events (If the events were to be come true, stop the ode34 function) ###
def event1(t, y):
    To, n_go, n_lo, r = y
    _, _, _, _, P_chmb, P = tank_system(t, To, n_go, n_lo, r)
    return P_chmb - P

event1.terminal = True
    # Event1 : Checks to see if Tank pressure and Chamber pressure are equal. When 0 stop

def event2(t, y):
    _, _, n_lo, _ = y
    return n_lo

event2.terminal = True
    # Event2 : Checks to see if the N20 Liquid mass in tank = 0

" SYSTEM OF EQUATIONS "
def odefun(t, y):
    # Set Initial conditions
    To, n_go, n_lo, r = y
    # Call ode45
    dT, dn_g, dn_l, dr, _ , _ = tank_system(t, To, n_go, n_lo, r)
    # Return changes in temperature, ullage gas, and liquid N20
    return [dT, dn_g, dn_l, dr]


def tank_system(t, To, n_go, n_lo, r):
    # Curve fitted combustion chamber pressure [Pa]:
    P_chmb = -2924.42 * t**6 + 46778.07 * t**5 - 285170.63 * t**4 + 813545.02 * t**3 - 1050701.53 * t**2 + 400465.85 * t + 1175466.2  # Test 1
    # Pe = 101325

    ### Given Constants ###
    n_He = 0.0                      # helium gas [kmol]
    Cd = 0.425                      # discharge coefficient: Test 1
    Ainj = 0.0001219352             # injector area [m**2]
    V = 0.0354                      # total tank volume [m**3]
    m_T = 6.4882                    # tank mass [kg]
    R = 8314.3                      # universal gas constant [J/(kmol*K)]
    MW2 = 44.013                    # molecular weight of N2O

    # Perry's Chemical Engineers' Handbook Property Equations
        # Temperature-dependent formulae which calculate specific heat capacities changes due to temperature
    G1 = 96.512                     # vapor pressure of N2O [Pa] coefficients
    G2 = -4045                      # valid for Temp range [182.3 K - 309.57 K]
    G3 = -12.277
    G4 = 2.886e-5
    G5 = 2
    Tc = 309.57                     # critical temperature of N2O [K]
    # heat of vaporization of N2O [J/kmol] coefficients
    J1 = 2.3215e7 
    J2 = 0.384                      # valid for Temp range [182.3 K - 309.57 K]
    J3 = 0
    J4 = 0
    # heat capacity of He at constant pressure [J/(kmol*K)] coefficients
    C1 = 0.2079e5                   
    C2 = 0                          # valid for Temp range [100 K - 1500 K]
    C3 = 0
    C4 = 0
    C5 = 0
    # heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients
    D1 = 0.2934e5
    D2 = 0.3236e5                   # valid for Temp range [100 K - 1500 K]
    D3 = 1.1238e3
    D4 = 0.2177e5
    D5 = 479.4
    # heat capacity of N2O liquid at constant pressure [J/(kmol*K)] coefficients
    E1 = 6.7556e4 
    E2 = 5.4373e1                   # valid for Temp range [182.3 K - 200 K]
    E3 = 0
    E4 = 0
    E5 = 0
    # molar specific volume of liquid N2O [m**3/kmol] coefficients
    Q1 = 2.781 
    Q2 = 0.27244
    Q3 = 309.57
    Q4 = 0.2882

    " BLOWDOWN CALCULATIONS "
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
    P = (n_He + n_go)*R*To / (V - n_lo*Vhat_l)                                                          # Calculating Oxidizer Tank Pressure
    a = m_T*Cp_T + n_He*CVhat_He + n_go*CVhat_g + n_lo*CVhat_l
    b = P*Vhat_l
    e = -delta_Hv + R*To
    f = -Cd*Ainj*np.sqrt(2/MW2)*np.sqrt((P-P_chmb)/Vhat_l)
    j = -Vhat_l*P_sat
    k = (V - n_lo*Vhat_l)*dP_sat
    m = R*To
    q = R*n_go

    Z=(-f*(-j*a + (q-k)*b)) / (a*(m+j) + (q-k)*(e-b))
    W=(-Z*(m*a + (q-k)*e)) / (-j*a + (q-k)*b)
    
    # Derivative Functions
    dT = (b*W+e*Z)/a
    dn_g = Z
    dn_l = W

    " Calculating Combustion Chamber "
    # Calculating oxidizer mass flow
    dm_ox = Cd * A * np.sqrt(2 * rho_ox * (P - P_chmb))

    ## Calculating Fuel Grain Regression
    A_port1 = 2 * np.pi * r * L
    dr = a * (dm_ox / A_port1)**n
    r = r + (dr * t)

    # Calculating fuel mass flow rate
    A_port = 2 * np.pi * r * L
    m_fuel = dr * t * A_port * rho_fuel

    # Calculating combustion chamber mass
    m_chmb = m_fuel + (dm_ox * t)
    rho_chmb = m_chmb / V_chmb

    # Returns
    # Change in Temperature [K], Mass of Nitrous Oxide Gas [kg], Mas of Nitrous Oxide Liquid [kg], Chamber Pressure [Pa], Tank Pressure [Pa]
    return dT, dn_g, dn_l, dr, P_chmb, P

### Solve the ODE system using solve_ivp (Initial Value Problem) ###
# Note: t_eval=np.linspace(0, 10, 500) will have 500 evenlyspaced time steps from 0 to 10 seconds
#       events=[event1, event2] passes the two defined events that will stop ode45 if conditions are met

sol = solve_ivp(odefun, [0, 10], y0, events=[event1, event2], t_eval=np.linspace(0, 10, 500))

# Extract the results
T = sol.y[0]                # Temperature over Time
n_g = sol.y[1]              # Mass of Nitrous Oxide Gas over Time
n_l = sol.y[2]              # Mass of Nitrous Oxide Liquid over Time
r = sol.y[3]                # Radius of Fuel Grain over Time

## Calculating pressure over time (sol.t)
# Creating vectors of zeros the length of time
P_chmb = np.zeros(len(sol.t))
P = np.zeros(len(sol.t))

# Iterating through the tank_system to back out P_chmb and P
for i, ti in enumerate(sol.t):
    _, _, _, _, P_chmb[i], P[i] = tank_system(ti, T[i], n_g[i], n_l[i], r[i])


"-------- PLOTTING ---------"
plt.figure()
plt.grid(True)
plt.plot(sol.t, T, 'r', linewidth=2)
plt.title('Temperature vs. Time')
plt.xlabel('Time [s]')
plt.ylabel('Temperature [K]')
plt.show()

plt.figure()
plt.grid(True)
plt.plot(sol.t, n_g * MW2, 'b', linewidth=2)
plt.plot(sol.t, n_l * MW2, 'g', linewidth=2)
plt.title('Mass of N20 vs. Time')
plt.xlabel('Time [s]')
plt.ylabel('Mass of N2O [kg]')
plt.legend(['Mass of N2O gas', 'Mass of N2O liquid'])
plt.show()

plt.figure()
plt.grid(True)
plt.plot(sol.t, P/6895, 'm', linewidth=2)
plt.plot(sol.t, P_chmb/6895, 'c', linewidth=2)
plt.title('Pressure vs. Time')
plt.xlabel('Time [s]')
plt.ylabel('Pressure [PSI]')
plt.legend(['tank pressure', 'chamber pressure'])
plt.show()
