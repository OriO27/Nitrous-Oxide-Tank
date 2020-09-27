"""
OOrtiz
V: 0.0

"""

from scipy.integrate import odeint
import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt


class Model(object):
    def __init__(self):
        #graviational acceleration
        global g
        g = 9.81 #m/s^2

        #tank outter radius
        self.ro = 0.1778/2 #m (7in) Diameter was input
        #tank inner radius
        self.ri = 0.14605/2 #m (5.75in) Diameter was input

        #density of wall material (aluminum for now) 
        self.rho_wall = 2710 #kg/m3
        #heat capacity of wall (al for now)
        self.C_wall = 897 #J/ kg k
        
        #length of the tank
        self.L_tank = 0.508 #m (~20in)
        #calculates tank volume
        self.V_tank = np.pi * (self.ro**2 - self.ri**2) * self.L_tank

        #percent of tank filled with liquid N2O
        self.fill_percent = 0 #% 

        #atmospheric temperature
        self.T_atm = 28 + 273.15 #K


        #Temperature of the Vapor part of wall
        self.T_wall_vap = self.T_atm #K | wall temperature at vapor region
        self.T_wall_liq = self.T_atm #K | wall temperature at liquid region
        self.T_lv_surf = 28 + 273.15 #K | inside temeprature at liquid vapor boundry
        self.T_liq = 28 + 273.15 #K | liquid temperature
        self.T_vap = 28 + 273.15 #K | vapor temperature
        
        
        N2O_mass = 15.0 #kg | Total mass of N2O
        fill_percent = 0 #% | percentage initially filled with liquid N2O
        t_stop = 30 * 60 #s (30min * 60sec) | total process time
        self.M_dot = N2O_mass / t_stop #kg/s |mass flow rate
        self.tspan = np.linspace(0, t_stop, 5000)
        self.Fill = False
        
        init_dim = self.Dim_Calcs((N2O_mass * fill_percent),PropsSI('D','T',self.T_wall_liq,'Q',0,'N2O'))        
        
        self.V_liq = init_dim[0]
        self.V_vap = init_dim[1]
        self.L_liq = init_dim[2]
        self.L_vap = init_dim[3]
        self.A_liq_in = init_dim[4]
        self.A_vap_in = init_dim[5]
        self.A_liq_out = init_dim[6]
        self.A_vap_out = init_dim[7]
        self.m_wall_liq = init_dim[8]
        self.m_wall_vap = init_dim[9]
 
        self.m_liq = N2O_mass * fill_percent
        self.rho_liq = PropsSI('D','T',self.T_wall_liq,'Q',0,'N2O')
        if self.Fill == False:
            self.M_dot = -self.M_dot #mass flow goes into the tank
            self.init_dim = self.Dim_Calcs((self.m_liq),PropsSI('D','T',self.T_wall_liq,'Q',0,'N2O'))
        
        self.z = [0, 0, 0, 0, 0] 
        self.x = odeint(self.ODE, self.z, self.tspan)

    def Dim_Calcs(self, m_liq, rho_liq):
        
        V_liq = m_liq / rho_liq #m^3 | volume of liquid in tank
        V_vap = self.V_tank - V_liq #m^3 | volume of vapor in tank
        L_liq = V_liq / (np.pi * self.ri**2) #m | Length of liquid region
        L_vap = V_vap / (np.pi * self.ri**2) #m | Length of vapor region
        if L_liq + L_vap > self.L_tank:
            print("Error in length calcs....")
        
        #Note::: Surface area does not include end caps; hollow tube assumed
        A_liq_in = 2 * np.pi * self.ri * L_liq #m^2 |inner surface area of tank liquid region
        A_vap_in = 2 * np.pi * self.ri * L_vap #m^2 |inner surface area of tank vapor region

        A_liq_out = 2 * np.pi * self.ro * L_liq #m^2 | outter surface area of tank liquid region
        A_vap_out = 2 * np.pi * self.ro * L_vap #m^2 | outter surface area of tank vapor region

        m_wall_liq = np.pi * (self.ro**2 - self.ri**2) * L_liq * self.rho_wall #kg | mass of wall in liqquid region
        m_wall_vap = np.pi * (self.ro**2 - self.ri**2) * L_vap * self.rho_wall #kg | mass of wall in liqquid region
        
        
        return V_liq, V_vap, L_liq, L_vap, A_liq_in, A_vap_in, A_liq_out, A_vap_out, m_wall_liq, m_wall_vap


    def thermo_air(self, T , P):
        #heat capacity of air at sea-level pressure and atmospheric temperature
        Cp_air = PropsSI('CP0MASS','P', P,'T', T ,'air') #J/kg/K
        #density of air
        rho_air = PropsSI('D','P', P ,'T', T,'air') #kg/m^3
        #air viscosity
        mu_air = PropsSI('V','P', P ,'T', T ,'air') #Pa s
        #air thermal conductivity
        k_air =  PropsSI('CONDUCTIVITY' , 'P', P, 'T', T, 'air') #W/m/K
        #Thermal Expansion Coefficient
        beta_air = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT' , 'P', P , 'T', T, 'air') #1/K
        return Cp_air, rho_air, mu_air, k_air, beta_air

    def Thermo_N2O_Liq(self, T, P):
        #PROPSSI Solving as if quality (x=0) for liquid properties... may not be accurate
        #Will also try having T as an input, or to calculate with P input and T input and average both values... 
        #Heat Capacity
        Cp_liq = PropsSI('CP0MASS','T', T,'Q', 0 ,'N2O') #J/kg/K
        #density 
        rho_liq = PropsSI('D','T', T ,'Q', 0,'N2O') #kg/m^3
        #Thermal Expansion Coefficient
        beta_liq = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT' , 'T', T , 'Q', 0, 'N2O') #1/K
        #Enthalpy
        enthalpy_liq = PropsSI('H', 'T', T, 'P', P, 'N2O') #J/kg
        #Saturation Enthalpy
        enthalpy_liq_sat = PropsSI('H', 'T', T, 'Q', 0, 'N2O') #J/kg 

        #The Following from... "Thermophysical properties of nitrous oxide" from IHS ESDU
        #Viscosity from EQ 4.9
        theta = ((309.59) -(5.24))/(T  -(5.24) )
        mu_liq = 0.0293423 * np.exp((1.6089) * (theta - 1) **(1/3) + 2.0439*(theta-1)**(4/3)) #mNs/m^2
        mu_liq = mu_liq/1000 #Ns/m^2
        # Thermal conductivity from EQ 4.11 
        # Note:: Literature suggests temp range ends at 10C, model below ends at T_crit... 
        # May Cause Instability in response.... 
        k_liq = 72.35 * (1 + 1.5 * (1 - (T/309.59))**(1/3) + (-3.5) * (1-(T/309.59))**(2/3) + 4.5 * (1- (T/309.59))) #mW/m/K
        k_liq = k_liq/1000

        return Cp_liq, rho_liq, mu_liq, k_liq, beta_liq, enthalpy_liq, enthalpy_liq_sat


    def Thermo_N2O_Vap(self, T, P):
        #PROPSSI Solving as if quality (x=0) for liquid properties... may not be accurate
        #Will also try having T as an input, or to calculate with P input and T input and average both values... 
        #Heat Capacity
        Cp_vap = PropsSI('CP0MASS','T', T,'Q', 1 ,'N2O') #J/kg/K
        #density 
        rho_vap = PropsSI('D','T', T ,'Q', 1,'N2O') #kg/m^3
        #Thermal Expansion Coefficient
        beta_vap = PropsSI('ISOBARIC_EXPANSION_COEFFICIENT' , 'T', T , 'Q', 1, 'N2O') #1/K
        #Enthalpy
        enthalpy_vap = PropsSI('H', 'T', T, 'P', P, 'N2O') #J/kg
        #Saturation Enthalpy
        enthalpy_vap_sat = PropsSI('H', 'T', T, 'Q', 1, 'N2O') #J/kg
        
        #The Following from... "Thermophysical properties of nitrous oxide" from IHS ESDU
        #Viscosity from EQ 4.9
        mu_vap = np.exp(3.3281 + (-1.18237) * (1/(T/309.59) - 1)**(1/3) + (-0.055155) * (1/(T/309.59) - 1)**(4/3)) #microNs/m^2
        mu_vap = mu_vap/ 1000000 #Ns/m^2
        #thermal conductivity
        k_vap =  np.exp(-7.0887 + (-0.276962) * (1-(T/309.59))**(-2/3) + 2.8872 * (1-(T/309.59))**(-1/3) + 16.6116 * \
                    (1-(T/309.59))**(1/3) + (-11.8221) * (1-(T/309.59))**(2/3)) #mW/m/K
        k_vap = k_vap/1000 #W/m/K
        return Cp_vap, rho_vap, mu_vap, k_vap, beta_vap, enthalpy_vap, enthalpy_vap_sat


        #z = [self.T_wall_vap, 0, 0,   self.T_wall_liq, self.T_lv_surf, self.T_liq, self.T_vap, self.V_liq, self.V_vap,\
        #    self.L_liq, self.L_vap, self.A_liq_in, self.A_vap_in, self.A_liq_out, self.A_vap_out

        #z = [0, 0, 0, 0, 0     ]

    def ODE(self, z, t):
        ###Need to rearrange initial conditions to match derivatives taken
        
        Q_liq_surf = z[0]        
        Q_surf_vap = z[1]        
        Q_wall_vap_out = z[2]  
        Q_wall_liq_out = z[3]        
        m_evap = z[4]

        print(Q_liq_surf)
        print(Q_surf_vap)
        print(Q_wall_vap_out)
        print(Q_wall_liq_out)
        print(m_evap)

        T_wall_vap = self.T_atm
        #Q_wall_vap_in = z[1]

        #Q_wall_vap_cond = z[3]
        #m_wall_vap_in = z[4]
        T_wall_liq = self.T_atm
        #Q_wall_liq_in = z[6]

        #Q_wall_liq_cond = z[8]
        #m_wall_liq_in = z[9]
        #m_vap = z[10]
        m_liq = self.m_liq



        #m_cond = z[15]
        T_vap = self.T_atm
        #U_vap = z[17]
        #rho_vap = 
        T_liq = self.T_atm
        #U_liq = z[20]
        #rho_liq = z[21]
        #Q_in_vap = z[23]
        #Q_in_liq = z[24]
        #P_vap = z[25]
        #P_liq = z[26]
        #m_out = z[27]

        dim = self.Dim_Calcs(m_liq, self.rho_liq)

        #V_liq = dim[0]
        #V_vap = dim[1]
        L_liq = dim[2]
        L_vap = dim[3]
        A_liq_in = dim[4]
        A_vap_in = dim[5]
        #A_liq_out = dim[6]
        #A_vap_out = dim[7]
        #m_wall_liq = dim[8]
        #m_wall_vap = dim[9]

        P = 101325 * 800 

        """
        #The film temperature between air and vapor wall region
        T_film_wall_vap = 0.5 * (z[0] + self.T_atm)
        #Finds air properties at film temp and atmospheric pressure
        Air = self.thermo_air(T_film_wall_vap, 101325)
        Cp_air = Air[0]
        rho_air = Air[1]
        mu_air = Air[2]
        k_air = Air[3]
        beta_air = Air[4]
        """
        T_film_liq_surf = 0.5 * (T_vap + T_liq)
        N2O_liq = self.Thermo_N2O_Liq(T_film_liq_surf, P)
        Cp_liq = N2O_liq[0]
        rho_liq = N2O_liq[1]
        mu_liq = N2O_liq[2]
        k_liq = N2O_liq[3]
        beta_liq = N2O_liq[4]
        enthalpy_liq = N2O_liq[5]
        enthalpy_liq_sat = N2O_liq[6]

        N2O_vap = self.Thermo_N2O_Vap(T_film_liq_surf, P)
        Cp_vap = N2O_vap[0]
        rho_vap = N2O_liq[1]
        mu_vap = N2O_vap[2]
        k_vap = N2O_vap[3]
        beta_vap = N2O_vap[4]
        #enthalpy_vap = N2O_vap[5]
        enthalpy_vap_sat = N2O_vap[6]


        N2O_latent_heat = enthalpy_vap_sat - enthalpy_liq_sat
        if L_liq <= 0: #a hacky fix to avoid dividing by 0 in a filling model
            L_liq = 0.00000000001 #m |very hacky approach.......

        EQ13 = (0.15 * ((Cp_liq * rho_liq**2 * g * beta_liq * abs(T_film_liq_surf - T_liq ) * L_liq**3 ) \
                / (mu_liq * k_liq ) )**(1/3) * (k_liq/L_liq)) * A_liq_in * (T_film_liq_surf - T_liq)

        EQ14 = (0.15 * ((Cp_vap * rho_vap**2 * g * beta_vap * abs(T_film_liq_surf - T_vap) * L_vap**3 )   \
                / (mu_vap * k_vap ) )**(1/3) * (k_liq/L_vap)) * A_vap_in * (T_film_liq_surf - T_liq)

        EQ2 = (0.021 * ((Cp_vap * rho_vap**2 * g * beta_vap * abs(T_wall_vap - T_vap) * L_vap**3  ) \
                / (mu_vap * k_vap ) )**(2/5) * (k_vap/L_vap)) * A_vap_in * (T_wall_vap - T_vap)

        EQ7 = (0.021 * ((Cp_liq * rho_liq**2 * g * beta_liq * abs(T_wall_liq - T_liq) * L_liq**3  ) \
                / (mu_liq * k_liq ) )**(2/5) * (k_liq/L_liq)) * A_liq_in * (T_wall_liq - T_liq)

        EQ12 = (EQ13 - EQ14) / (N2O_latent_heat + (enthalpy_liq_sat - enthalpy_liq) )

        return EQ13, EQ14, EQ2, EQ7, EQ12



Model()




"""
        dz[0] = (dz[1] - dz[2] + dz[3] + dz[4] * self.C_wall * (self.T_wall_vap - z[0]) )/ \
                self.m_wall_vap * self.C_wall 

        dz[1] = (0.59 * ((Cp_air * (rho_air )**2 * g * beta_air * abs(z(1) - self.T_atm) * z(30)**3  )\
                / (mu_air * k_air ) )**(0.25) * (k_air/z(30))) * A_wall_vap_out * (z(1) - self.T_atm)

        dz[2] = (0.021 * ((Cp_vap * (z(19) )^2 * g * biot_vap * abs(z(1) - z(17)) * z(30)^3  )\
                / (mu_vap * k_vap )   )^(2/5) * (k_vap/z(30))) * A_wall_vap_inside * (z(1) - z(17))



        self.m_wall_vap += dz[4]
        
"""

