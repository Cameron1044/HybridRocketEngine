#Converted IdealLiquidModel.m to Python

import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
m_loaded = 19.32933  # N2O mass initially loaded into tank [kg]: Test 1
Ti = 286.5  # initial temperature [K]: Test 1
V = 0.0354  # total tank volume [m**3]
R = 8314.3  # universal gas constant [J/(kmol*K)]
MW2 = 44.013  # molecular weight of N2O
G1 = 96.512  # vapor pressure of N2O [Pa] coefficients
G2 = -4045  # valid for Temp range [182.3 K - 309.57 K]
G3 = -12.277
G4 = 2.886e-5
G5 = 2
Q1 = 2.781  # molar specific volume of liquid N2O [m**3/kmol] coefficients
Q2 = 0.27244
Q3 = 309.57
Q4 = 0.2882

n_to = m_loaded / MW2  # initial total N2O in tank [kmol]
Vhat_li = Q2**(1 + (1 - Ti / Q3)**Q4) / Q1  # molar volume of liquid N2O [m**3/kmol]
To = Ti  # initial temperature [K]
P_sato = np.exp(G1 + G2 / To + G3 * np.log(To) + G4 * To**G5)  # initial vapor pressure of N20 [Pa]
n_go = P_sato * (V - Vhat_li * n_to) / (-P_sato * Vhat_li + R * To)  # initial N2O gas [kmol]
n_lo = (n_to * R * To - P_sato * V) / (-P_sato * Vhat_li + R * To)  # initial N2O liquid [kmol]

# Initial conditions
y0 = [Ti, n_go, n_lo]

def event1(t, y):
    To, n_go, n_lo = y
    _, _, _, Pe, P = tank_system(t, To, n_go, n_lo)
    return Pe - P

event1.terminal = True

def event2(t, y):
    _, _, n_lo = y
    return n_lo

event2.terminal = True

def odefun(t, y):
    To, n_go, n_lo = y
    dT, dn_g, dn_l, _, _ = tank_system(t, To, n_go, n_lo)
    return [dT, dn_g, dn_l]


def tank_system(t, To, n_go, n_lo):
    Pe = -2924.42 * t**6 + 46778.07 * t**5 - 285170.63 * t**4 + 813545.02 * t**3 - 1050701.53 * t**2 + 400465.85 * t + 1175466.2  # Test 1

    # Given constants
    n_He = 0 # helium gas [kmol]
    Cd = 0.425 # discharge coefficient: Test 1
    Ainj = 0.0001219352 # injector area [m**2]
    V = 0.0354 # total tank volume [m**3]
    m_T = 6.4882 # tank mass [kg]
    R = 8314.3 # universal gas constant [J/(kmol*K)]
    MW2 = 44.013 # molecular weight of N2O
    G1 = 96.512 # vapor pressure of N2O [Pa] coefficients
    G2 = -4045 # valid for Temp range [182.3 K - 309.57 K]
    G3 = -12.277
    G4 = 2.886e-5
    G5 = 2
    Tc = 309.57 # critical temperature of N2O [K]
    J1 = 2.3215e7 # heat of vaporization of N2O [J/kmol] coefficients
    J2 = 0.384 # valid for Temp range [182.3 K - 309.57 K]
    J3 = 0
    J4 = 0
    C1 = 0.2079e5 # heat capacity of He at constant pressure [J/(kmol*K)] coefficients
    C2 = 0 # valid for Temp range [100 K - 1500 K]
    C3 = 0
    C4 = 0
    C5 = 0
    D1 = 0.2934e5 # heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients
    D2 = 0.3236e5 # valid for Temp range [100 K - 1500 K]
    D3 = 1.1238e3
    D4 = 0.2177e5
    D5 = 479.4
    E1 = 6.7556e4 # heat capacity of N2O liquid at constant pressure [J/(kmol*K)] coefficients
    E2 = 5.4373e1 # valid for Temp range [182.3 K - 200 K]
    E3 = 0
    E4 = 0
    E5 = 0
    Q1 = 2.781 # molar specific volume of liquid N2O [m**3/kmol] coefficients
    Q2 = 0.27244
    Q3 = 309.57
    Q4 = 0.2882

    Vhat_l = Q2**(1+(1-To/Q3)**Q4)/Q1
    CVhat_He = C1 + C2*To + C3*To**2 + C4*To**3 + C5*To**4 - R
    CVhat_g = D1 + D2*((D3/To)/np.sinh(D3/To))**2 + D4*((D5/To)/np.cosh(D5/To))**2 - R
    CVhat_l = E1 + E2*To + E3*To**2 + E4*To**3 + E5*To**4
    Tr = To/Tc
    delta_Hv = J1*(1 - Tr) ** (J2 + J3*Tr + J4*Tr**2)
    P_sat = np.exp(G1 + G2/To + G3*np.log(To) + G4*To**G5)
    dP_sat = (-G2/(To**2) + G3/To + G4*G5*To**(G5-1)) * np.exp(G1 + G2/To + G3*np.log(To) + G4*To**G5)
    Cp_T = (4.8 + 0.00322*To)*155.239
    P = (n_He + n_go)*R*To / (V - n_lo*Vhat_l)

    a = m_T*Cp_T + n_He*CVhat_He + n_go*CVhat_g + n_lo*CVhat_l
    b = P*Vhat_l
    e = -delta_Hv + R*To
    f = -Cd*Ainj*np.sqrt(2/MW2)*np.sqrt((P-Pe)/Vhat_l)
    j = -Vhat_l*P_sat
    k = (V - n_lo*Vhat_l)*dP_sat
    m = R*To
    q = R*n_go

    Z=(-f*(-j*a + (q-k)*b)) / (a*(m+j) + (q-k)*(e-b))
    W=(-Z*(m*a + (q-k)*e)) / (-j*a + (q-k)*b)
    
    dT = (b*W+e*Z)/a
    dn_g = Z
    dn_l = W

    return dT, dn_g, dn_l, Pe, P


# Solve the ODE system using solve_ivp
sol = solve_ivp(odefun, [0, 10], y0, events=[event1, event2], t_eval=np.linspace(0, 10, 500))

# Extract the results
T = sol.y[0]
n_g = sol.y[1]
n_l = sol.y[2]

Pe = np.zeros(len(sol.t))
P = np.zeros(len(sol.t))
for i, ti in enumerate(sol.t):
    _, _, _, Pe[i], P[i] = tank_system(ti, T[i], n_g[i], n_l[i])

# Plotting code
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
plt.plot(sol.t, P, 'm', linewidth=2)
plt.plot(sol.t, Pe, 'c', linewidth=2)
plt.title('Pressure vs. Time')
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
plt.legend(['tank pressure', 'chamber pressure'])
plt.show()
