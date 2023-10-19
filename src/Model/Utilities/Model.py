import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from src.Model.Utilities.Blowdown import BlowdownModel
from src.Model.Utilities.Combustion import CombustionModel

class Model():
    def __init__(self, initialInputs):
        self.iterations = 500
        self.tspan = [0, 10]
        self.y0 = [0, 0, 0, 0, 0]
        self.blowdown = BlowdownModel()
        self.combustion = CombustionModel()

    def event1(self):
        value = 0

    def event2(self):
        value = 0

    def model(self, mo, mf, Po, Pc, Vu, r):
        dmo_dt, dPo_dt, dVu_dt = self.blowdown.blowdownModel(Po, Vu, r, mf)
        dr_dt, dmf_dt, dPc_dt = self.combustion.combustionModel(dmo_dt, Pc, r)
        return [dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt]

    def ODE45(self):
        sol = solve_ivp(self.model, self.tspan, self.y0, events=[self.event1, self.event2], t_eval=np.linspace(self.tspan[0], self.tspan[1], self.iterations))

