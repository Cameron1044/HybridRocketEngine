# File for Nozzle Model
# being given an array of each time step for champer pressure, mass flow rate of fuel and oxidizer, tank pressure, temperature of tank and chamber

# https://www.sciencedirect.com/science/article/pii/S2405844018374164  justification for using alpha =15 and beta =45

import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Read in necessary files
#   Read in chamber pressure array from some kind of ZK Output File
p_chmb_arr = range(1047, 500, 50)
print(p_chmb_arr)

# Constants
gamma = 1.249 #Ratio of Specific Heats
R = 188.9 # [J/kg*K] gas constant for nitrous oxide
Pa = 14.98 # [psi] Ambient Pressure (Assumed in Boulder, Colorado)
mp = 0 #mass flow of propellant 
mo = 0 #mass flow of oxidizer
C_str = 0 #characteristic velocity
Ptc = 0 #Chamber Total Pressure and can be total pressure Pt for nozzle stations
mdot = mo + mp #mass flow of oxidizer and propellant
alpha = 15 # [Degrees] - typical half cone angle for majority of nozzle designs
beta = 45 #degrees -typical converging angle for nozzle diffuser
Pe = 0 #exit pressure assumed to be ambient pressure for perfectly expanded nozzle
Tc = 0; #chamber temperature same as total temperature Tt in nozzle


# Alt function, calculates ideal area ratio for each timestep

# Preallocate arrays
Me_arr = []
AreaRatios = []

for i in range (0,len(p_chmb_arr)): # inefficient for loop method but idc research more later lol

    # Calculate total pressure over exit pressure
    P0oPa = p_chmb_arr(i)/Pa

    # Backsolve for exit mach number
    Me_arr.append(np.sqrt(((P0oPa)**((gamma - 1) / gamma) - 1 )/( (gamma - 1) / 2 )))

    # Use area-mach relation to find Ae/A* (ideally = Ae/At), with many intermediate values to simplify code
    intval1 = (gamma+1)/(gamma-1)
    intval2 = 2/(gamma+1)
    intval3 = 1 - (((gamma-1)/2) * Me_arr(i)**2)
    AreaRatios.append(np.sqrt(1/(Me_arr(i)**2) * (intval2 * intval3)**intval1))

    # Function to define nozzle geometry:
def Nozzle_Characteristics(mdot,C_str,Ptc,Pe,Dc, AeoAt,alpha,beta,R):
    At = (mdot*C_str)/Ptc # Calculates optimal throat area given mass flow rate, characteristic velcoity, and chamber pressure
    
    Me = np.sqrt((2/(gamma-1))*((Ptc/Pe)**((gamma-1)/gamma)-1))
    Te = Tc/(1+((gamma-1)/2)*(Me**2)) #exit temperature 
    Pthr = Ptc*(2/(gamma+1))**(gamma/(gamma-1)) #Pressure at throat of nozzle assuming Mach = 1 at throat
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
 # Currently missing Pe, the exit pressure at the end of the nozzle. May very well be a problem like the nozzle problems from aero :/ ***Pe is the same as Pa*** @jacob
    return T #Thrust is returned, to be used to develop thrust profile

At, rt, re, Length_cone, Ae, AeoAt, Dt, De, Me, Te, Pthr, Length_diff, Ve = Nozzle_Characteristics(mdot, C_str, Ptc, Pe, Tc, gamma, alpha, beta, R)

# Displaying the results
print("Throat Area (At):", At)
print("Throat Radius (rt):", rt)
print("Exit Radius (rt):", re)
print("Cone Length in mm (Length_cone):", Length_cone)  # Note: Length_cone is repeated twice in the return statement.
print("Exit Area (Ae):", Ae)
print("Area Ratio (AeoAt):", AeoAt)
print("Throat Diameter (Dt):", Dt)
print("Exit Diameter (De):", De)
print("Mach Number at Exit (Me):", Me)
print("Exit Temperature (Te):", Te)
print("Pressure at Throat (Pthr):", Pthr)
print("Length of Diffuser (Length_diff):", Length_diff)
print("Exit Velocity (Ve):", Ve)
# test
# Plot Thrust Profile? 