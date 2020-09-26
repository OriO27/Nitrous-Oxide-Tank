"""
OOrtiz
V: 0.0

"""

import scipy.integrate as spi
import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP

class Model():
    def __init__(self):
        #graviational acceleration
        self.g = 9.81 #m/s^2

        #tank outter radius
        self.ro = 0.127 #m (5in)

        #tank inner radius
        self.ri = 0.127 #m (5in)

        #density of wall material (aluminum for now) 
        self.rho_wall = 2,710 #kg/m3
        #heat capacity of wall (al for now)
        self.C_wall = 897 #J/ kg k
        
        #length of the tank
        self.L_tank = 0.508 #m (~20in)
        
        #calculates tank volume
        self.V_tank = np.pi * (ro**2 - ri**2) * self.L_tank

        #percent of tank filled with liquid N2O
        self.fill_percent = 0 #% 

        #atmospheric temperature
        self.T_atm = 28 + 273.15 #K


        #Temperature of the Vapor part of wall
        
        self.T_wall_vap = 30 + 273.15 #K


    def thermo_air(self):
        #heat capacity of air at sea-level pressure and atmospheric temperature
        self.Cp_air = PropsSI('CP0MASS','P',101325,'T', self.T_atm ,'air') #J/kg/K
        #density of air
        self.rho_air = PropsSI('D','P',101325,'T',self.T_atm,'air') #kg/m^3
        #air viscosity
        self.mu_air = PropsSI('V','P',101325,'T',(300),'air') #Pa s
        #air thermal conductivity
        self.k_air =  PropsSI('CONDUCTIVITY' , 'P', 101325, 'T', self.T_atm, 'air') #W/m/K
    
    def Initial_conditions(self, z):
        #assumes initial temperature is atm temperature 
        z[0] = self.T_atm


        #mass of the vapor region wall
        z[4] = self.rho_wall * ###VOLUME OF VAPOR REGION WALL

    def ODE(self, t, dz):

        dz[0] = (dz[1] - dz[2] + dz[3] + dz[4] * self.C_wall * (self.T_wall_vap - z[0]) )/ \
                self.m_wall_vap * self.C_wall 


        dz[1] = ()

            (0.59 * ((self.Cp_air * (self.rho_air )^2 * g * biot_air * abs(z(1) - self.T_atm) * z(30)^3  )...
                / (self.mu_air * k_air )   )^(0.25)...
                * (k_air/z(30))) * A_wall_vap_out * (z(1) - self.T_atm)


        self.m_wall_vap += dz[4]
        
