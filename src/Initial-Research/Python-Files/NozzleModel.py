# File for Nozzle Model
# being given an array of each time step for champer pressure, mass flow rate of fuel and oxidizer, tank pressure, temperature of tank and chamber

# https://www.sciencedirect.com/science/article/pii/S2405844018374164  justification for using alpha =15 and beta =45

import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
gamma = 1.249 #change to correct gamma
R = 8.3144 #J/mol*k
Pa = 0 # Ambient Pressure in Boulder, Colorado
mp = 0 #mass flow of propellant 
mo = 0 #mass flow of oxidizer
C_str = 0 #characteristic velocity
Pc = 0 #Chamber Pressure and can be total pressure Pt for nozzle stations
mdot = mo + mp #mass flow of oxidizer and propellant
alpha = 15 # [Degrees] - typical half cone angle for majority of nozzle designs
beta = 45 #degrees -typical converging angle for nozzle diffuser
Pe = 0 #exit pressure assumed to be ambient pressure for perfectly expanded nozzle
Tc = 0; #chamber temperature same as total temperature Tt in nozzle

# Function to define nozzle geometry:
def Nozzle_Characteristics(mdot,C_str,Pc,Pe,Dc, AeoAt,alpha,beta):
    At = (mdot*C_str)/Pc # Calculates optimal throat area given mass flow rate, characteristic velcoity, and chamber pressure
    
    Me = np.sqrt((2/(gamma-1))*((Pc/Pe)**((gamma-1)/gamma)-1))
    Te = Tc/(1+((gamma-1)/2)*(Me**2)) #exit temperature 
    Pthr = Pc*(2/(gamma+1))**(gamma/(gamma-1)) #Pressure at throat of nozzle assuming Mach = 1 at throat
    AeoAt = (1/Me)*((2/(gamma+1))*(1+((gamma-1)/2)*(Me**2)))**((gamma+1)/(2*(gamma-1)))
    Ae = AeoAt * At # Calculates exit area w/ throat area and area ratio (given constants)
    rt = np.sqrt(At/np.pi) # Determines throat radius based on throat area
    re = np.sqrt(Ae/np.pi) # Determines exit radius based on exit area
    Dt = 2 * rt # Throat Diameter
    De = 2 * re # Exit Diameter
    Length_cone = (0.5 * (De - Dt)) / (np.tan(alpha)) # Calculate cone length parallel to x axis
    Length_diff = (0.5 * (Dc - Dt)) / (np.tan(beta)) #calculate lenght of diffuser parallel to x axis
    Ve = Me*np.sqrt(gamma*R*Te)
    return At, rt, Length_cone, Ae,Dt,De,Length_cone,Me,Te,Pthr, Length_diff,Ve

# Function to predict thrust performance w/ nozzle geometry
def CalculateThrust(mdot,Ve, Ae, Pe):
    T = (mdot * Ve) + ((Pe - Pa) * Ae)
 # Currently missing Pe, the exit pressure at the end of the nozzle. May very well be a problem like the nozzle problems from aero :/
    return T #Thrust is returned, to be used to develop thrust profile

# test
# Plot Thrust Profile? 