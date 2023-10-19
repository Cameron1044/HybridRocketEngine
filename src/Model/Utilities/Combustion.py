import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class CombustionModel():
    def __init__(self):
        self.a = 0
        self.n = 0
        self.L = 0
        self.Ru = 0
        self.T_chmb = 0
        self.M_chmb = 0
        self.gamma = 0
        self.rho_f = 0
        self.A_t = 0
        
    def fuelGrainState(self, r, Pc):
        # Take necessary constant valus from init dunder function
        Ru = self.Ru
        L = self.L
        M_chmb = self.M_chmb
        T_chmb = self.T_chmb

        # Calculate the current state of the combustion chamber and fuel grain
        R_prime = Ru / M_chmb
        A_p = np.pi * r**2
        A_b = 2 * np.pi * r * L
        V_chmb = A_p * L
        rho_chmb = Pc / (R_prime * T_chmb)
        return A_p, A_b, V_chmb, rho_chmb, R_prime

    def fuelRegression(self, dmo_dt, A_p):
        # Take necessary constant values from init dunder function
        a = self.a
        n = self.n

        # Define Derivative for Fuel Grain port radius from Marxman Fuel Regression
        dr_dt = a * (dmo_dt / A_p)**n
        return dr_dt

    def mfDeriv(self, dr_dt, A_b):
        # Take necessary constant values from init dunder function
        rho_f = self.rho_f
        
        # Define Derivative for Fuel Mass Flow from Marxman Fuel Regression
        dmf_dt = -A_b * rho_f * dr_dt
        return dmf_dt

    def PChmbDeriv(self, A_b, V_chmb, rho_chmb, R_prime, dr_dt, dmo_dt, Pc):
        # Take necessary constant values from init dunder function
        rho_f = self.rho_f
        A_t = self.A_t
        gam = self.gamma
        M_chmb = self.M_chmb
        T_chmb = self.T_chmb

        # Calculating Fuel Mass relations from gas generation and fuel regression
        FuelMassRelation = A_b * dr_dt * (rho_f - rho_chmb)
        # Calculating mass flow through a Nozzle
        NozzleMass = Pc * A_t * np.sqrt(gam / (R_prime * T_chmb)) * (2 / (gam + 1))^((gam + 1)/(2 * (gam - 1))); 
        # Define Derivative for Oxidizer Mass Flow from Bernoulli's Equation
        dPc_dt = (R_prime * T_chmb/V_chmb) * (FuelMassRelation - dmo_dt - NozzleMass);   
        return dPc_dt

    def combustionModel(self, dmo_dt, Pc, r):
        A_p, A_b, V_chmb, rho_chmb, R_prime = self.fuelGrainState(r, Pc)
        dr_dt = self.fuelRegression(dmo_dt, A_p)
        dmf_dt = self.mfDeriv(dr_dt, A_b)
        dPc_dt = self.PChmbDeriv(A_b, V_chmb, rho_chmb, R_prime, dr_dt, dmo_dt, Pc)
        return dr_dt, dmf_dt, dPc_dt