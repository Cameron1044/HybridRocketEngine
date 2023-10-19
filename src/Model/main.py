import numpy as np
import math as m
import matplotlib.pyplot as plt
from tabulate import tabulate
import pandas as pd

# from Utilities.Model import Model
from Utilities.utilities import units

initialInputs = { # Dictionary of initial inputs
    # Injector
    "A_inj": units(9.4e-6, 'm^2'),
    "C_d": units(0.4, 'unitless'),
    # Fuel Properties
    "rho_fuel": units(919, 'kg/m^3'),
    "M": units(23.147, 'g/mol'),
    "T": units(3015, 'K'),
    "gamma": units(1.2705, 'unitless'),
    "n": units(0.681, 'unitless'),
    "a": units(2.85E-5, 'a'),
    # Fuel Grain
    "L_fuel": units(12, 'in'),
    "OD_fuel": units(1.688, 'in'),
    "ID_fuel": units(1.35, 'in'),
    # Oxidizer Tank
    "V_tank": units(3, 'L'),
    "P_tank": units(3000, 'psi'),
    "ullage_fraction": units(0.20, 'unitless'),
    # Nozzle
    "d_t": units(0.5, 'in'),
    # Ambient Conditions
    "P_amb": units(102675.3, 'Pa'),

    ##### CONSTANTS #####
    "rho_ox": units(1226, 'kg/m^3'),
    "Ru": units(8.3143, 'J/(mol*K)'),
}

print(initialInputs)
# model = Model(initialInputs)