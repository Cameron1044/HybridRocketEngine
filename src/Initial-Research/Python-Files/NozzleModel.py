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
AeoAt = 3.31
alpha = 15 # [Degrees] - typical half cone angle for majority of nozzle designs

# Function to define nozzle geometry:
def Nozzle_Characteristics():
    At = ((mp+mo)*C_str)/Pc # Calculates optimal throat area given mass flow rate, characteristic velcoity, and chamber pressure
    rt = np.sqrt(At/np.pi) # Determines throat radius based on throat area
    Ae = AeoAt * At # Calculates exit area w/ throat area and area ratio
    re = np.sqrt(Ae/np.pi) # Determines exit radius based on exit area

    # Calculate nozzle length
    Dt = 2 * rt # Throat Diameter
    De = 2 * re # Exit Diameter
    Length_cone = (1/2 * (De - Dt)) / (np.tan(alpha))

    return At, rt, Length_cone

