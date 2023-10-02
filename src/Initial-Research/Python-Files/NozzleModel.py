# File for Nozzle Model
#To Do: set up function to take inputs for formulas to calculate outputs

# being given an array of each time step for champer pressure, mass flow rate of fuel and oxidizer, tank pressure, temperature of tank and chamber

import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

mp = 0
C_str = 0
Pc = 0
mo = 0

def Nozzle_Characteristics():
    At = ((mp+mo)*C_str)/Pc


    return

