import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class ChemicalProperties():
    def __init__(self, gas="He", mass_loaded=0, tank_volume=0, initial_temp=0):
        self.To = initial_temp # initial temperature [K]
        self.mass_loaded = mass_loaded # N2O mass initially loaded into tank [kg]
        self.V = tank_volume # total tank volume [m**3]
        self.MW_N2O = 44.013 # molecular weight of N2O
        self.Tc = 309.57 # critical temperature of N2O [K]
        self.R = 8314.3 # universal gas constant [J/(kmol*K)]

        C = self.ullageGas(gas) # heat capacity of gas at constant pressure [J/(kmol*K)] coefficients
        self.C1 = C[0]
        self.C2 = C[1]
        self.C3 = C[2]
        self.C4 = C[3]
        self.C5 = C[4]

        self.G1 = 96.512  # vapor pressure of N2O [Pa] coefficients
        self.G2 = -4045  # valid for Temp range [182.3 K - 309.57 K]
        self.G3 = -12.277
        self.G4 = 2.886e-5
        self.G5 = 2

        self.Q1 = 2.781  # molar specific volume of liquid N2O [m**3/kmol] coefficients
        self.Q2 = 0.27244
        self.Q3 = 309.57
        self.Q4 = 0.2882

        self.J1 = 2.3215e7 # heat of vaporization of N2O [J/kmol] coefficients
        self.J2 = 0.384 # valid for Temp range [182.3 K - 309.57 K]
        self.J3 = 0
        self.J4 = 0

        self.D1 = 0.2934e5 # heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients
        self.D2 = 0.3236e5 # valid for Temp range [100 K - 1500 K]
        self.D3 = 1.1238e3
        self.D4 = 0.2177e5
        self.D5 = 479.4

        self.E1 = 6.7556e4 # heat capacity of N2O liquid at constant pressure [J/(kmol*K)] coefficients
        self.E2 = 5.4373e1 # valid for Temp range [182.3 K - 200 K]
        self.E3 = 0
        self.E4 = 0
        self.E5 = 0

    def ullageGas(self, gas):
        C1_N2 = 0.2911e5 # heat capacity of N2 at constant pressure [J/(kmol*K)] coefficients
        C2_N2 = 0.0861e5
        C3_N2 = 1.7016e3
        C4_N2 = 0.0010e5
        C5_N2 = 909.79

        C1 = 0.2079e5 # heat capacity of He AND argon at constant pressure [J/(kmol*K)] coefficients
        C2 = 0 # valid for Temp range [100 K - 1500 K]
        C3 = 0
        C4 = 0
        C5 = 0
        gasDict = {
            "N2O": [C1_N2, C2_N2, C3_N2, C4_N2, C5_N2],
            "He": [C1, C2, C3, C4, C5],
            "Ar": [C1, C2, C3, C4, C5]
        }
        return gasDict[gas]
    
    def initialMoles(self):
        To = self.To  # initial temperature [K]
        n_to = self.mass_loaded / self.MW_N2O  # initial total N2O in tank [kmol]
        Vhat_li = self.Q2**(1 + (1 - To / self.Q3)**self.Q4) / self.Q1  # molar volume of liquid N2O [m**3/kmol]
        P_sato = np.exp(self.G1 + self.G2 / To + self.G3 * np.log(To) + self.G4 * To**self.G5)  # initial vapor pressure of N20 [Pa]
        n_go = P_sato * (self.V - Vhat_li * n_to) / (-P_sato * Vhat_li + self.R * To)  # initial N2O gas [kmol]
        n_lo = (n_to * self.R * To - P_sato * self.V) / (-P_sato * Vhat_li + self.R * To)  # initial N2O liquid [kmol]
        return n_go, n_lo
    
    def chemicalStates(self, T):
        Tr = T/self.Tc
        Vhat_l = self.Q2**(1+(1-T/self.Q3)**self.Q4)/self.Q1
        CVhat_He = self.C1 + self.C2*T + self.C3*T**2 + self.C4*T**3 + self.C5*T**4 - self.R
        CVhat_g = self.D1 + self.D2*((self.D3/T)/np.sinh(self.D3/T))**2 + self.D4*((self.D5/T)/np.cosh(self.D5/T))**2 - self.R
        CVhat_l = self.E1 + self.E2*T + self.E3*T**2 + self.E4*T**3 + self.E5*T**4
        delta_Hv = self.J1*(1 - Tr) ** (self.J2 + self.J3*Tr + self.J4*Tr**2)
        P_sat = np.exp(self.G1 + self.G2/T + self.G3*np.log(T) + self.G4*T**self.G5)
        dP_sat = (-self.G2/(T**2) + self.G3/T + self.G4*self.G5*T**(self.G5-1)) * np.exp(self.G1 + self.G2/T + self.G3*np.log(T) + self.G4*T**self.G5)
        Cp_T = (4.8 + 0.00322*T)*155.239
        return Vhat_l, CVhat_He, CVhat_g, CVhat_l, delta_Hv, P_sat, dP_sat, Cp_T
