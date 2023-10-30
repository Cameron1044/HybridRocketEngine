import numpy as np
import math as m
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from .ChemicalProperties import ChemicalProperties

class ZKmodel(ChemicalProperties):
    def __init__(self, inputs):
        self.A_inj = inputs["A_inj"]                                # Area of the Injector
        self.C_d = inputs["C_d"]                                    # Coefficient of Discharge of Injector
        self.P_tank = inputs["P_tank"]                                # Pressure of the Oxidizer Tank
        self.V_tank = inputs["V_tank"]                              # Volume of the Oxidizer Tank
        self.T_tank = inputs["T_tank"] # initial temperature [K]
        self.m_T = inputs["m_T"] # total mass of tank [kg]
        self.m_N2O = inputs["m_N2O"] # mass of loaded N2O [kg]
        super().__init__(T_tank=self.T_tank, V_tank=self.V_tank, P_tank=self.P_tank, m_N2O=self.m_N2O)
    
    def ZKModel(self, To, n_go, n_lo, P_chmb):

        Vhat_l = self.molarVolumeN2O(To)
        CVhat_Ar = self.specificHeatArgon(To)
        CVhat_g = self.specificHeatN2OGas(To)
        CVhat_l = self.specificHeatN2OLiquid(To)
        delta_Hv = self.heatVaporizationN2O(To)
        P_sat = self.vaporPressureN2O(To)
        dP_sat = self.vaporPressureN2ODerivative(To)
        Cp_T = self.specificHeatAluminum(To)
        
        Po = self.tankPressure(n_lo, n_go, self.n_Ar, To, self.V_tank)

        a = self.m_T*Cp_T + self.n_Ar*CVhat_Ar + n_go*CVhat_g + n_lo*CVhat_l
        b = Po*Vhat_l
        e = -delta_Hv + self.R*To
        f = -self.C_d*self.A_inj*np.sqrt(2/self.MW_N2O)*np.sqrt((Po-P_chmb)/Vhat_l)
        j = -Vhat_l*P_sat
        k = (self.V_tank - n_lo*Vhat_l)*dP_sat
        m = self.R*To
        q = self.R*n_go

        Z=(-f*(-j*a + (q-k)*b)) / (a*(m+j) + (q-k)*(e-b))
        W=(-Z*(m*a + (q-k)*e)) / (-j*a + (q-k)*b)
        
        # Derivative Functions
        dT = (b*W+e*Z)/a
        dn_g = Z
        dn_l = W

        return dT, dn_g, dn_l, Po
