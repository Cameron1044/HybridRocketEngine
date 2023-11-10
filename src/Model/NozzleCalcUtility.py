import numpy as np
import math as m
import pandas as pd
import matplotlib.pyplot as plt


## Inputs/Givens

# Pressures
Pc =  362 # Pressure of the Combustion Chamber [psi] 
Pc_arr = np.linspace(Pc, 100)
P0 = Pc # Since flow is essentially stagnant in combustion chamber, we consider this the total pressure of the flow
Pamb = 14.93 # Ambient pressure at launch site (Fort Collins, CO)
P0oPe = P0/Pamb

# Geometry
Dt = 1.2 # Diameter of the Throat [in]
At = np.pi/4 * Dt**2 # Area of the Throat (Calculated) [in^2]

# Misc
gamma = 1.205 # Ratio of speicifc heats

## Calculate all the things 

# Calculate exit mach number from pressure ratio
Me = np.sqrt( ( (P0oPe) ** ((gamma-1)/gamma) - 1) / ( (gamma - 1) / 2) )
print("Mach at Nozzle exit is", Me)

# Armed with exit Mach number, calcualte area ratio and exit area 
Arat = np.sqrt((1/Me**2)*((2/(gamma + 1))*(1 + (gamma-1)*Me**2/2))**((gamma+1)/(gamma-1))) 
Ae = Arat * At
print("Correct Nozzle Area Ratio is", Ae)
De = np.sqrt(4*Ae/np.pi)
print("Correct Nozzle Exit Diameter is", De, "inches")

# Ratio of P0 to Pe is PRESERVED regardless of how P0 changes, just a quirk of having choked flow
Pe_arr = []
for pressures in Pc_arr:
    Pe_arr.append((P0oPe)**(-1) * pressures)

PeoPamb = []
for exitpressures in Pe_arr:
    PeoPamb.append(exitpressures/Pamb)

for pressureRatios in PeoPamb:
    if pressureRatios < 0.4:
        print("Approximate Chamber Pressure for Summerfield Flow Separation Condition is ", Pc_arr[PeoPamb.index(pressureRatios)])
        break

# Calculate Thrust and Thrust Coefficient for Given Nozzle Design

# Preallocate empty arrays
Cf = []
F = []

# Fill arrays with a for loop
for index, Pe in enumerate(Pe_arr):

    Cf.append(np.sqrt( (2*(gamma)**2/gamma - 1) * (2/(gamma+1))**((gamma+1)/(gamma-1)) * (1 - (Pe_arr[index]/Pc_arr[index])**((gamma - 1)/gamma))))

    F.append(At * Pc_arr[index] * Cf[index])

# Plot force as a function of time / chamber pressure
plt.figure()
plt.plot(F, Pc_arr)
plt.show()

