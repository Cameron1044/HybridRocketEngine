import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Initial conditions
Cd = 0.175 #Coefficient of discharge, orifice | dimensionless
n_holes = 4 #Number of holes, injector | quantity
phi = 0.0492 #Diameter of injection holes | inches
rho_ox = 0.0440 #Density of oxidizer | lb/in^3
p_ox = 2300 #Pressure of oxidizer tank | psi
p_chmb = 0 #Pressure of Combustion Chamber, pre-fire | psi
V_ox = 15.8213 #Total volume of oxidizer tank | in^3
L = 4 #Fuel port length | inches
r_port0 = 0.25 #Initial fuel port radius | inches
a = 0.7 #Burn rate coefficient | scalar
n = 0.8 #Pressure exponent | scalar
d_t = 0.23 #Throat diameter | inches
gamma = 1.249
cd_throat = 0.2
Ru = 8.314
T = 3059 #1
M = 28.811

# Constants
A = n_holes*np.pi*(phi/2)**2 #Area of injection hole | in^2
print(A)
mdot_ox0 = Cd*A*np.sqrt(2*rho_ox*(p_ox - p_chmb)*32.2*12) #Initial mass flow rate of oxidizer | lb/s
V_ox_liquid = V_ox*(80/100) #Initial volume of liquid oxidizer | in^3
V_ox_gas0 = V_ox - V_ox_liquid #Initial volume of gas pressurant | in^3
P1V1 = p_ox*V_ox_gas0 #Initial pressure * volume in oxidizer tank

rho_fuel = 0.0421 #Density of fuel | lb/in^3
A_t = np.pi*(d_t/2)**2 #Cross-sectional area of nozzle throat | in^2

tf=100
tstep=0.031
i_f=tf/tstep

# Initial
p_ox = 2300
p_chmb = 0
V_ox_gas = V_ox_gas0
Vox_liquid = 0
A_port = 2 * np.pi * r_port0 * L
r = r_port0

dr_dt = 0
dm_ox_dt = mdot_ox0
dV_dt = 0

t_arr = []
p_chmb_arr = []
thrust_choked_arr = []

def tank_system(t, To, n_go, n_lo, Pe, Ainj, Cd, n_He, V):
    ### Given Constants ###
    # n_He = 0.0                      # helium gas [kmol]
    # Cd = 0.425                      # discharge coefficient: Test 1
    # Ainj = 0.0001219352             # injector area [m**2]
    # V = 0.0354                      # total tank volume [m**3]
    m_T = 6.4882                    # tank mass [kg]

    R = 8314.3                      # universal gas constant [J/(kmol*K)]
    MW2 = 44.013                    # molecular weight of N2O
    # Perry's Chemical Engineers' Handbook Property Equations
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

    ### BLOWDOWN CALCULATIONS ###
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
    Cp_T = (4.8 + 0.00322*To)*155.239                                                                   # specific heat of tank, Aluminum [J/(kg*K)]
    
    ## Simplified expression definitions for solution
    P = (n_He + n_go)*R*To / (V - n_lo*Vhat_l)
    # print(P, Pe, n_He, n_go, n_lo)
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
    
    # Derivative Functions
    dT = (b*W+e*Z)/a
    dn_g = Z
    dn_l = W

    return dn_g, dn_l, dT, P

def initial(m_loaded, V, Ti):
    R = 8314.3                              # universal gas constant [J/(kmol*K)]
    MW2 = 44.013                            # molecular weight of N2O
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

    ### Calculating Initial Conditions ###
    n_to = m_loaded / MW2                                                   # initial total N2O in tank [kmol]
    Vhat_li = Q2**(1 + (1 - Ti / Q3)**Q4) / Q1                              # molar volume of liquid N2O [m**3/kmol]
    To = Ti                                                                 # initial temperature [K]
    P_sato = np.exp(G1 + G2 / To + G3 * np.log(To) + G4 * To**G5)           # initial vapor pressure of N20 [Pa]
    print(P_sato/6895, Vhat_li, R, To, n_to)
    n_go = P_sato * (V - Vhat_li * n_to) / (-P_sato * Vhat_li + R * To)     # initial N2O gas [kmol]
    n_lo = (n_to * R * To - P_sato * V) / (-P_sato * Vhat_li + R * To)      # initial N2O liquid [kmol]
    return n_go, n_lo

Ainj = A/1550
# print(A)
Cd = 0.175
n_He = 0.000125
V = V_ox/61020
To = 294.261
m_loaded = 0.4/2.205
n_go, n_lo = initial(m_loaded, V, To)
# print(n_go, n_lo)
p_chmb = 101325/6895
n_lo_arr = []
n_go_arr = []

for i in range(0,int(i_f)):
    if p_ox < p_chmb or n_lo <= 0:
        break

    dn_g, dn_l, dT, p_ox = tank_system(i, To, n_go, n_lo, p_chmb*6895, Ainj, Cd, n_He, V)
    To = To + dT*tstep
    n_go = n_go + dn_g*tstep
    n_lo = n_lo + dn_l*tstep
    n_lo_arr.append(n_lo)
    n_go_arr.append(n_go)

    MW2 = 44.013
    dm_ox_dt = n_lo*MW2*2.205
    p_ox = p_ox/6895

    # dm_ox_dt = Cd * A * np.sqrt(2 * rho_ox * (p_ox - p_chmb) * 32.2 * 12)
    # dV_dt = dm_ox_dt / rho_ox * tstep
    # V_ox_gas = V_ox_gas + dV_dt
    dr_dt = a * (dm_ox_dt / A_port)**n #regression

    r = r + dr_dt*tstep

    A_port = 2 * np.pi * r * L
    m_fuel = dr_dt*tstep*A_port*rho_fuel
    p_chmb = (3.28084*(dm_ox_dt + m_fuel/tstep))*(840.8786447/0.588693023)/(32.2*A_t)

    t_arr.append(i*tstep)
    p_chmb_arr.append(p_chmb)

    m_ox = dm_ox_dt*tstep
    m_chmb = m_fuel + m_ox
    V_chmb = np.pi*(1.25/2)**2*0.5 + np.pi*(1.25/2)**2*1 + np.pi*(r)**2*L
    rho_chmb = m_chmb / V_chmb
    mdot_exit_choked = cd_throat * A_t * np.sqrt(32.2*12*gamma*rho_chmb*p_chmb*(2/(gamma+1))**((gamma+1)/(gamma-1)))
    v_exit = 3.28084*np.sqrt(((1000*2*gamma)/(gamma-1))*((Ru*T)/M))
    thrust_choked = mdot_exit_choked * v_exit / 32.2

    thrust_choked_arr.append(thrust_choked)

    if np.isnan(dm_ox_dt) or np.isnan(dr_dt) or np.isnan(dV_dt) or np.isnan(p_chmb) or np.isnan(thrust_choked):
        break

    print(dm_ox_dt, r, p_chmb, p_ox, thrust_choked)

plt.figure()
plt.grid(True)
plt.plot(t_arr, p_chmb_arr)
plt.xlabel('Time (s)')
plt.ylabel('P_chmb (psi)')
plt.show()

plt.figure()
plt.grid(True)
plt.plot(t_arr, thrust_choked_arr)
plt.xlabel('Time (s)')
plt.ylabel('Thrust (psi)')
plt.show()

plt.figure()
plt.grid(True)
plt.plot(t_arr, n_go_arr, 'b', linewidth=2)
plt.plot(t_arr, n_lo_arr, 'g', linewidth=2)
plt.title('Mass of N20 vs. Time')
plt.xlabel('Time [s]')
plt.ylabel('Mass of N2O [kg]')
plt.legend(['Mass of N2O gas', 'Mass of N2O liquid'])
plt.show() # chris kongy boy