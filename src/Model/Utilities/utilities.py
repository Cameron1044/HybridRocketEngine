import pandas as pd
import matplotlib.pyplot as plt

def units(value, conversion, n=1):
    ## To take in a value and the current unit it is in to change it into metric
    ## INPUTS: value - the value we want to convert
    ##         unit  - the unit that will be changed to metric
    ##         n     - Fuel Regression constant used to calculate a
    ##
    ## OUTPUTS: The Value will have it's units changed into metric / SI units
    ##          If the unit is already in metric, nothing will happen
    conversionDict = {
        ## For Units that are English
        # Length --> m
        "in": 1/39.37,
        "ft": 1/3.281,
        # Surfrace Area --> m2
        "in^2": 1/1550,
        "ft^2": 1/10.764,
        # Volume --> m3
        "in^3": 1/61020,
        "ft^3": 1/35.315,
        # Pressure --> Pa
        "psi": 6895,
        # Mass --> kg
        "lbm": 1/2.205,
        # Force --> N
        "lbf": 4.448,
        # Density --> kg/m^3
        "lbm/ft^3": 16.01846337396,

        ## For Units that are Metric
        # Length --> m
        "m": 1,
        "cm": 1/100,
        "mm": 1/1000,
        # Surfrace Area --> m2
        "m^2": 1,
        "cm^2": 1/10000,
        "mm^2": (1e-6),
        # Volume --> m3
        "m^3": 1,
        "cm^3": (1e-6),
        "mm^3": (1e-9),
        "L": 1/1000,
        # Pressure --> Pa
        "Pa": 1,
        "kPa": 1000,
        # Mass --> kg
        "kg": 1,
        "g": 1/1000,
        # Force --> N
        "N": 1,
        "kN": 1000,
        # Density --> kg/m^3
        "kg/m^3": 1,
        # Molecular Weight --> kg/mol
        "g/mol": 1/1000,
        # Temperature --> K
        "K": 1,
        # Gas Constant --> J/(mol*K)
        "J/(mol*K)": 1,

        ## Unitless
        "unitless": 1,

        ## For a Burn coefficient
        "a": (0.0254**(1 + 2*(n)))*(0.453592**(-n))
    }
    return value * conversionDict[conversion]