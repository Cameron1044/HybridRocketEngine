import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class CombustionModel():
    def __init__(self):

        self.rho_f = 0 # Density of Fuel Material
        self.a = 0 # Marxman Constant #1 (Fuel Dependant)
        self.n = 0 # Marxman Constant #2 (Fuel Dependant)
        self.R = 0 # Specific Gas Constant (Fuel Dependent)
        self.Tc = 0 # Chamber Temperature [K]
        self.At = 0 # Throat Area [m^2]
    
    def BurnRate(self, dmo_dt, A_port): # Variables that change throughout burn time
        # Take necessary constant values from init dunder function
        a = self.a
        n= self.n

        # Define Burn Rate from Marxman's Equation
        dr_dt = -a*(dmo_dt/A_port)**(n)

        return dr_dt
    
    def FuelMassFlow(self, A_b, dr_dt): # Variables that change throughout burn time
        # Take necessary constant values from init dunder function
        rho_f = self.rho_f

        # Define Fuel Mass Flow Derivative
        dmf_dt = -A_b*rho_f*dr_dt
        return dmf_dt

    def ChamberPressureDeriv(self, dr_dt, Ab, Pc, Vc, dmo_dt):

        # Take necessary constant values from init dunder function
        rho_f = self.rho_f
        R = self.R
        Tc = self.Tc
        At = self.At
        gamma = self.gamma
        # Intermediate calculation of "Chamber Density", Fuel Accumulation Rate, and Nozzle Mass Flow Rate
        rho_c = Pc/(R*Tc)
        fuel = Ab*dr_dt*(rho_f-rho_c)
        dm_nozzle_dt = Pc*At*np.sqrt(gamma/(R*Tc))*(2/(gamma+1))**((gamma+1)/(2*(gamma-1))); 

        # Calculate Derivative of Chamber Pressure 
        dPc_dt = R*(Tc/Vc) * (fuel - dmo_dt - dm_nozzle_dt); 

        return dPc_dt
    
    def combustionModel(self, dmo_dt, Pc, A_port, Ab, Vc):

        # Calculate All Combustion Derivatives
        dr_dt =  self.BurnRate(self,dmo_dt, A_port)
        dmf_dt = self.FuelMassFlow(self, Ab, dr_dt)
        dPc_dt = self.ChamberPressureDeriv(self, dr_dt, Ab, Pc, Vc, dmo_dt)

        return dmf_dt, dr_dt, dPc_dt