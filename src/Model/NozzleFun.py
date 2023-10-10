import numpy as np
from scipy.optimize import fsolve


"----- Aerodynamic Functions -----"
# The following are isentropic and nozzle functions used in ASEN 3111 from Professor Jensens Aerodynamic code base converted into Python
# They will be used in the main nozzle script function
def A_A_star_M_solve(gamma, AoAstar):
    ## Python function which will solve the subsonic and supersonic mach numbers for a given gamma and A over A star area ratio
    def star(x): # Area-Mach Equation
        return ((1/x**2) * ((2/(gamma+1)) * (1 + ((gamma-1)/2) * x**2))**((gamma+1)/(gamma-1))) - AoAstar**2

    Msup = fsolve(star, 2, xtol=1e-8, maxfev=1000)
    Msub = fsolve(star, 0.1, xtol=1e-8, maxfev=1000)
    # Similar to MATLAB fsolve however the options in MATLAB need to be explicit in Python for functionality

    return Msup, Msub

def A_A_star_solve(gamma, MachNumber):
    ## Python function which will solve A over A star area ratio given a gamma and Mach Number
        # Note: Does not matter if MachNumber is supersonic or subsonic. Will return same AoAstar
    AoAstar = np.sqrt((1/MachNumber**2) * ((2 / (gamma + 1)) * (1 + (gamma - 1) * MachNumber**2 / 2))**((gamma + 1) / (gamma - 1)))
    return AoAstar

def isentropic(gamma, M):
    ## Python function to calculate isentropic flow relationships between total intrinsic properties and static intrinsic properites
    g = gamma
    p0op = (1 + (gamma - 1) / 2 * M**2) ** (gamma / (gamma - 1))
    t0ot = 1 + (gamma - 1) / 2 * M**2

    return p0op, t0ot

def isentropicFindM(gamma, p0op):
    ## Python function to solve a Mach Number from a ratio of total pressure to local pressure and gamma
    Msub = np.sqrt(((p0op)**((gamma -1) / gamma) - 1) * 2 / (gamma - 1))
    return Msub

