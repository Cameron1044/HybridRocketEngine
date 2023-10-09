import numpy as np
import math as m

"----- Author's Notes -----"
# In reality, the oxidizer tank cannot be treated as a Ideal Gas Law and Incompressible flow problem.
# In fact, a self-pressurizing system is not in phase equillibrium because there are irreversibilities in heat transer and work
# However, this script will hopefully give ball parks on the initial conditions of the oxidizer tank

## PURPOSE:
# Based on Adrian's volume of 3-L
    # OUTPUTS: Mass Breakdown, Moles of Liquid Oxidizer and Moles of Ullage gas
"----- Inital Conditions -----"
"Oxidizer Tank - Liquid Nitroux Oxide"
## Liquid Nitroux Oxide
R = 8.3143                              # universal gas constant [J/(mol*K)]
MW2 = 44.013                            # molecular weight of N2O [g / mol]
rho_ox = 1226                           # Liquid Nitroux Oxide at 273K [kg/m³]
MWGN = 28.014                           # Molecular weight of N2 [g/mol]
MWGA = 40                               # Molecular weight of Ar [g/mol]

## Tank
Ti_tank = 298                           # Initial temperature of the Oxidizer tank | K
V_tank = 0.003                          # Oxidizer Tank Volume | [m^3] (3L tank)
P_ox = 2500 * 6895                      # Initial pressume of the oxidizer tank | [Psi --> Pa]

"----- Calculations -----"
## Splitting Volume of Oxidizer Tank
V_ox_liquid = V_tank*(80/100)           # Initial volume of liquid oxidizer | [L or 1000 cm^3 --> m^3]
    # Note - assuming 80% of the volume is liquid nitrous oxide
V_gas  = V_tank*(15/100)                # Initial volume of Ullage Gas | [m^3]
    # Note - assuming 18% of the remaining volume is for ullage gas
    # Thus 2% of the volume is for Phase Equillibirum between Nitrous Gas liquid and gas phases

## Incompressible Liquid Nitrous Oxide
m_loaded = V_ox_liquid * rho_ox         # Initial mass of liquid oxidizer | kg
n_liquid = m_loaded / MW2 * 1000
print("\n")
print("Number of Moles of Liquid Nitroux Oxide (mol): ", n_liquid)
print("Weight of Liquid Oxidizer (kg): ", m_loaded)

## Ideal Gas Law of Ullage Tank
# PV = nRT solving for number of moles
n_gas = (P_ox * V_gas) / (R * Ti_tank)
m_gas = n_gas * MWGN / 1000
print("Number of Moles of Ullage gas (mol): ", n_gas)
print("Weight of Ullage gas (kg): ", m_gas)

## Total System
m_system = m_loaded + m_gas
print("Total System Weight (kg): ", m_system)
print("Total System weight (lb): ", m_system * 2.205, "\n")
