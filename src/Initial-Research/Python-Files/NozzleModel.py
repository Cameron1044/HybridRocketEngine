# File for Nozzle Model
# being given an array of each time step for champer pressure, mass flow rate of fuel and oxidizer, tank pressure, temperature of tank and chamber

import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
gc = 32.2 # [ft/s^2] - Gravitational constant ONLY USE IF NOT WORKING IN SI UNITS
Pa = 0 # Ambient Pressure
mp = 0 #mass flow of propellant 
mo = 0 #mass flow of oxidizer
C_str = 0
Pc = 0 #Chamber Pressure
mdot = mo + mp #mass flow of oxidizer and propellant
AeoAt = 3.31 #Area ratio
alpha = 15 # [Degrees] - typical half cone angle for majority of nozzle designs

# Function to define nozzle geometry:
def Nozzle_Characteristics(mdot,C_str,Pc,AeoAt,alpha):
    At = (mdot*C_str)/Pc # Calculates optimal throat area given mass flow rate, characteristic velcoity, and chamber pressure
    rt = np.sqrt(At/np.pi) # Determines throat radius based on throat area
    Ae = AeoAt * At # Calculates exit area w/ throat area and area ratio (given constants)
    re = np.sqrt(Ae/np.pi) # Determines exit radius based on exit area

    # Calculate nozzle length
    Dt = 2 * rt # Throat Diameter
    De = 2 * re # Exit Diameter
    Length_cone = (1/2 * (De - Dt)) / (np.tan(alpha)) # Calculate cone length

    return At, rt, Length_cone, Ae

# Function to predict thrust performance w/ nozzle geometry
def CalculateThrust(mdot,Ve, Ae, Pe):
    T = (mdot * Ve / gc) + ((Pe - Pa) * Ae)
 # Currently missing Pe, the exit pressure at the end of the nozzle. May very well be a problem like the nozzle problems from aero :/
    return T #Thrust is returned