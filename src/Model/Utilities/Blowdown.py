import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class BlowdownModel():
    def __init__(self):
        self.A_i = 0 # Injector Area (m^2)
        self.C_d = 0 # Injector Coefficient of Discharge
        self.rho_ox = 0 # Density of Oxidizer (N20)
        self.Po_i = 0 # Initial Chamber Pressure 
        self.V_u_i = 0 # Initial Ullage Volume 


    def OxidizerMassDeriv(self, Po, P_c): # Changing inputs are tank pressure and chamber pressure
        
        # Take necessary constant values from init dunder function
        A_i = self.A_i
        C_d = self.C_d
        rho_ox = self.rho_ox

        # Define Derivative for Oxidizer Mass Flow from Bernoulli's Equation
        dmo_dt = -A_i*C_d*np.sqrt(2*rho_ox*(Po - P_c))
        return dmo_dt

    def VolDeriv(self, dmo_dt):
        # Take necessary constant values from init dunder function
        rho_ox = self.rho_ox

        # Define Derivative of Ullage Volume from Oxidizer Mass Flow Rate
        dVu_dt = -dmo_dt/rho_ox

        return dVu_dt

    def TankPressureDeriv(self, Vu):

        # Take necessary constant values from init dunder function
        Po_i = self.Po_i
        V_u_i = self.V_u_i

        # Define Derivative of Tank Pressure from Ideal Gas & constant Temperature assumption
        dP0_dt = -((Po_i * V_u_i)/Vu**2)

        return dP0_dt
    
    def BlowdownModel(self, Po, Pc, Vu, dmo_dt):

        # Calculate all blowdown derivatives
        dmo_dt = self.OxidizerMassDeriv(self,Po,Pc)
        dVu_dt = self.VolDeriv(self,dmo_dt)
        dP0_dt = self.TankPressureDeriv(self, Vu)

        return dmo_dt, dVu_dt, dP0_dt
