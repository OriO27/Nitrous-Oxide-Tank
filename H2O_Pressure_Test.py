"""
OOrtiz 
Pressure Vessel Evacuation
V:: 0.0
"""

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from scipy import *
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

class N2OModel():

    def __init__(self):
        self.Cd = .9 #coefficient of friction at injector site
        self.A_i = .003 #[m^2] total cross sectional area of the injector holes
        self.T0 = 306 #[K] initial temperature 
        self.P0 = 650000 #[Pa] initial pressure
        self.t = 0.01 #[s] timestep ##will change to an integral approach in later iterations
        self.V_tank = 0.07413327496 #[m^3]
        self.m_tot = 30 #[kg] intial total tank mass
        vapPercent = 0.15
        liqPercent = 0.85
        m_dot_0 = self.startCond(vapPercent, liqPercent)
        self.propChange(self.T0, self.m_tot, m_dot_0)      



    def startCond(self, vapPercent, liqPercent ):
        Rho_l = PropsSI('D', 'T', self.T0, 'Q', 1 , 'N2O')
        Rho_v = PropsSI('D', 'T', self.T0, 'Q', 0 , 'N2O')        
        P1 = PropsSI('P', 'T', self.T0, 'Q', 1, 'N2O')
        self.m_l_0 = self.m_tot * liqPercent #[kg] initial liquid mass percent   
        m_dot = self.mFlux(Rho_l, P1)
        return m_dot

    def mFlux(self, Rho1, P1): #function calculates current mass flux, Pressure input should be at saturation
        m_spi = self.A_i * (self.Cd * np.sqrt(2*Rho1 * P1))
        return m_spi




    def propChange(self, T_curr, m_tot, m_dot):
        """
        From CoolProp.org:::
        The latent heat of vaporization can be calculated using the difference between the vapor and liquid enthalpies at the same point on the saturation curve.
        """

        m_curr = m_tot - m_dot * self.t
        m_l_old = self.m_tot - m_dot * self.t

        Rho_l = PropsSI('D', 'T', T_curr, 'Q', 1 , 'N2O')
        Rho_v = PropsSI('D', 'T', T_curr, 'Q', 0 , 'N2O') 

        m_l_new = (self.V_tank - m_curr/Rho_v)

        h_v = PropsSI('H', 'T', T_curr, 'Q', 0, 'N2O')
        h_l = PropsSI('H', 'T', T_curr, 'Q', 1, 'N2O')
        h_vapor = h_v - h_l
        dQ = m_dot * h_vapor
        C = PropsSI('C', 'T', T_curr, 'Q', 1, 'N2O')
        dT = -dQ / (m_l_old * C)
        print(dT)
        print(C)







N2OModel()