# Import necessary libraries
import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd

# Import custom modules
from .Blowdown import BlowdownModel
from .Combustion import CombustionModel
from .utilities import ToEnglish, ToMetric

class Model():
    ##### ##### ##### ##### ##### ##### ##### ##### #####
    ## This class outlines the entire Rocket Model     ##
    ##### ##### ##### ##### ##### ##### ##### ##### #####

    def __init__(self, initialInputs, iterations=1000, tspan=[0, 20]):
        # Purpose:  Initiates the class
        # Inputs:   initialInputs - The input dictionary from main.py
        #           iterations    - Number of data points across the time span
        #           tspan         - Array containing the start and finish time for interation
        #
        # Outputs:  self          - Input Constants necessary to perform the Calculations

        # Comes from Main.py input dictionary and command call
        self.initialInputs = initialInputs 
        self.iterations = iterations
        self.tspan = tspan
        
        #### Calculate initial conditions based on provided inputs ####
        # Mass of Liquid Nitrous Oxide
        mo_i = initialInputs["V_tank"] * (1 - initialInputs["ullage_fraction"]) * initialInputs['rho_ox']
        # Mass of Fuel Grain
        mf_i = initialInputs['rho_fuel'] * np.pi * initialInputs['L_fuel'] * (initialInputs['OD_fuel']**2 - initialInputs['ID_fuel']**2) / 4
        # Initial Oxidizer Tank Pressure
        Po_i = initialInputs["P_tank"]
        # Inital Chamber Pressure
        Pc_i = initialInputs["P_amb"]
        # Initial Volume of the Oxidizer Tank Ullage Gas
        Vu_i = initialInputs["V_tank"] * initialInputs["ullage_fraction"]
        # Initial Port Radius of the Fuel Grain
        r_i = initialInputs['ID_fuel']/2
        self.mt_i = mo_i + mf_i

        ## Loading into ODE45 initial variable array
        self.y0 = [mo_i, mf_i, Po_i, Pc_i, Vu_i, r_i, 0]
        
        # Initialize sub-models for blowdown and combustion
        self.blowdown = BlowdownModel(initialInputs)
        self.combustion = CombustionModel(initialInputs)

    #### TERMINATION EVENTS ####
    ## Purpose: Events which stop the integration and iteration if they become true

    # Event to check if port radius exceeds half of the fuel outer diameter
    def event1(self, t, y):
        r = y[5]
        fuel_OD = self.initialInputs['OD_fuel']
        return r - fuel_OD/2

    # Event to check if oxidizer mass reaches zero
    def event2(self, t, y):
        m_o = y[0]
        return m_o

    # Event to check if chamber pressure exceeds oxidizer tank pressure
    def event3(self, t, y):
        P_c = y[3]
        P_o = y[2]
        return P_c - P_o

    def model(self, mo, mf, Po, Pc, Vu, r, I):
        # Purpose:  Initiates the Model which wraps around the different Hybrid rocket Models
        
        # Calls the Blowdown Model
        dmo_dt, dPo_dt, dVu_dt = self.blowdown.blowdownModel(Po, Pc, Vu)

        # Calls the Combustion Chamber Model
        dr_dt, dmf_dt, dPc_dt, T = self.combustion.combustionModel(dmo_dt, Pc, r)

        return dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt, T

    # Wrapper function for the ODE solver
    def ODEfun(self, t, y):
        mo, mf, Po, Pc, Vu, r, I = y
        dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt, T = self.model(mo, mf, Po, Pc, Vu, r, I)
        return [dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt, T]

    ## Function to solve the ODEs and return results in a dataframe
    # This is where the CSV files will be made and how 
    def ODE45(self):
        # Define columns for the results dataframe
        data_columns = [
            "time",
            "thrust",
            "impulse",
            "Pc",
            "Pox",
            "dmox",
            "dmf",
            "OF",
            "r",
            "dm_total_dt",
            "cstar",
        ]

        # Initialize dataframe
        self.df = pd.DataFrame(columns=data_columns)
        
        # Solve the ODEs
        sol = solve_ivp(self.ODEfun, self.tspan, self.y0, t_eval=np.linspace(self.tspan[0], self.tspan[1], self.iterations), events=[self.event1, self.event2, self.event3])

        # Extract results from the solution
        t = sol.t
        mo = sol.y[0]
        mf = sol.y[1]
        Po = sol.y[2]
        Pc = sol.y[3]
        Vu = sol.y[4]
        r = sol.y[5]
        I = sol.y[6]

        # Populate the dataframe with results
        for i, ti in enumerate(sol.t):
            dmo_dt, dmf_dt, dPo_dt, dPc_dt, dVu_dt, dr_dt, T = self.model(mo[i], mf[i], Po[i], Pc[i], Vu[i], r[i], I[i])
            new_data = {
                "time": t[i],
                "thrust": ToEnglish(T, 'N'),
                "impulse": ToEnglish(I[i], 'N'),
                "Pc": ToEnglish(Pc[i], 'Pa'),
                "Pox": ToEnglish(Po[i], 'Pa'),
                "dmox": ToEnglish(-dmo_dt, 'kg'),
                "dmf": ToEnglish(-dmf_dt, 'kg'),
                "OF": ToEnglish(dmo_dt/dmf_dt, 'unitless'),
                "r": ToEnglish(r[i], 'm'),
                "dm_total_dt": ToEnglish(-(dmo_dt+dmf_dt), 'kg'),
                "cstar": ToEnglish(T, 'N')/ToEnglish(-(dmo_dt+dmf_dt), 'kg')
            }
            self.df.loc[len(self.df)] = new_data

        return self.df

# Set the terminal attribute for the events
Model.event1.terminal = True
Model.event2.terminal = True
Model.event3.terminal = True
