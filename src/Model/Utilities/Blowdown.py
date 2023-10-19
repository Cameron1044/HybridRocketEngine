import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class BlowdownModel():
    def __init__(self, inputs):
        self.A_inj = inputs["A_inj"]
        self.C_d = inputs["C_d"]
        self.rho_ox = inputs["rho_ox"]
        self.Po_i = inputs["P_tank"]
        self.V_u_i = inputs["V_tank"] * inputs["ullage_fraction"]

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

    def PressureDeriv(self, Vu):

        # Take necessary constant values from init dunder function
        Po_i = self.Po_i
        V_u_i = self.V_u_i

        # Define Derivative of Tank Pressure from Ideal Gas & constant Temperature assumption
        dP0_dt = -((Po_i * V_u_i)/Vu**2)

        return dP0_dt
    
    def BlowdownModel(self, Po, Pc, Vu, dmo_dt):
        dmo_dt = self.OxidizerMassDeriv(Po, Pc)
        dVu_dt = self.VolDeriv(dmo_dt)
        dP0_dt = self.PressureDeriv(Vu)
        return dmo_dt, dVu_dt, dP0_dt
