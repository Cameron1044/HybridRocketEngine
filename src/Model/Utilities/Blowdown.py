import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class BlowdownModel():
    def __init__(self, input1=0):
        self.input1 = input1
        self.A_i = 0
        self.C_d = 0
        self.rho_ox = 0
        self.P_0_i = 0
        self.V_u_i = 0

    def IdealGas(self):
        number = self.input1**2
        return dP_dt

    def OxidizerMassDeriv(self, P_0, P_c): # Changing inputs are tank pressure and chamber pressure
        
        # Take necessary constant values from init dunder function
        A_i = self.A_i
        C_d = self.C_d
        rho_ox = self.rho_ox

        # Define Derivative for Oxidizer Mass Flow from Bernoulli's Equation
        dmo_dt = -A_i*C_d*np.sqrt(2*rho_ox*(P_0 - P_c))
        return dmo_dt

    def VolDeriv(self, dmo_dt):
        # Take necessary constant values from init dunder function
        rho_ox = self.rho_ox

        # Define Derivative of Ullage Volume from Oxidizer Mass Flow Rate
        dVu_dt = -dmo_dt/rho_ox

        return dVu_dt

    def PressureDeriv(self, V):

        # Take necessary constant values from init dunder function
        P_0_i = self.P_0_i
        V_u_i = self.V_u_i

        # Define Derivative of Tank Pressure from Ideal Gas & constant Temperature assumption
        dP0_dt = -((P_0_i * V_u_i)/V**2)

        return dP0_dt
    
    def BlowdownModel(self):
        dmo_dt = self.OxidizerMassDeriv()
        dVu_dt = self.VolDeriv()
        dP0_dt = self.PressureDeriv()
        return dmo_dt, dVu_dt, dP0_dt
