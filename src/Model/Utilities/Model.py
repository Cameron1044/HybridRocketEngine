import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd

from .Blowdown import BlowdownModel
from .Combustion import CombustionModel
from .utilities import ToEnglish

class Model():
    def __init__(self, initialInputs, iterations=500, tspan=[0, 10]):
        self.initialInputs = initialInputs
        self.iterations = iterations
        self.tspan = tspan
        mo_i = initialInputs["V_tank"] * (1 - initialInputs["ullage_fraction"]) * initialInputs['rho_ox']
        mf_i = initialInputs['rho_fuel'] * np.pi * (initialInputs['ID_fuel']**2 - initialInputs['ID_fuel']**2) / 4
        Po_i = initialInputs["P_tank"]
        Pc_i = initialInputs["P_amb"]
        Vu_i = initialInputs["V_tank"] * initialInputs["ullage_fraction"]
        r_i = initialInputs['ID_fuel']/2
        self.y0 = [mo_i, mf_i, Po_i, Pc_i, Vu_i, r_i, 0]
        self.blowdown = BlowdownModel(initialInputs)
        self.combustion = CombustionModel(initialInputs)

    def event1(self, t, y):
        r = y[5]
        fuel_OD = self.initialInputs['OD_fuel']
        return r - fuel_OD/2

    def event2(self, t, y):
        m_o = y[0]
        return m_o

    def event3(self, t, y):
        P_c = y[3]
        P_o = y[2]
        return P_c - P_o

    def model(self, mo, mf, Po, Pc, Vu, r, I):
        dmo_dt, dPo_dt, dVu_dt = self.blowdown.blowdownModel(Po, Pc, Vu)
        dr_dt, dmf_dt, dPc_dt, T = self.combustion.combustionModel(dmo_dt, Pc, r)
        return dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt, T

    def ODEfun(self, t, y):
        mo, mf, Po, Pc, Vu, r, I = y
        dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt, T = self.model(mo, mf, Po, Pc, Vu, r, I)
        return [dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt, T]

    def ODE45(self):
        data_columns = [ # Dataframe columns
            "time",
            "thrust",
            "impulse",
            "Pc",
            "Pox",
            "dmox",
            "dmf",
            "OF",
            "r",
        ]

        self.df = pd.DataFrame(columns=data_columns) # Create dataframe
        sol = solve_ivp(self.ODEfun, self.tspan, self.y0, t_eval=np.linspace(self.tspan[0], self.tspan[1], self.iterations), events=[self.event1, self.event2, self.event3])

        t = sol.t
        mo = sol.y[0]
        mf = sol.y[1]
        Po = sol.y[2]
        Pc = sol.y[3]
        Vu = sol.y[4]
        r = sol.y[5]
        I = sol.y[6]

        for i, ti in enumerate(sol.t):
            dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt, T = self.model(mo[i], mf[i], Po[i], Pc[i], Vu[i], r[i], I[i])
            new_data = { # Dictionary of new data
                "time": t[i],
                "thrust": ToEnglish(T, 'N'),
                "impulse": ToEnglish(I[i], 'N'),
                "Pc": ToEnglish(Pc[i], 'Pa'),
                "Pox": ToEnglish(Po[i], 'Pa'),
                "dmox": ToEnglish(-dmo_dt, 'kg'),
                "dmf": ToEnglish(-dmf_dt, 'kg'),
                "OF": ToEnglish(dmo_dt/dmf_dt, 'unitless'),
                "r": ToEnglish(r[i], 'm'),
            }
            self.df.loc[len(self.df)] = new_data # Append new data to dataframe

        return self.df

Model.event1.terminal = True
Model.event2.terminal = True
Model.event3.terminal = True