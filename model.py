"""
OOrtiz
V: 1.1

Modeling Equations from 
"Review and Evalution of Models for 
Self-Pressurizing Propellant Tank Dynamics"
Zimmerman, J.E., Waxman, B.S., Cantwell, B.J.,
Zilliac, G.G. 49th Joint Propulsion Conference, 2013  

"""

from scipy.integrate import odeint
import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import os 

class Model(object):
    def __init__(self):
        self.clear = os.system('cls')
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
        self.k_wall = 150 #W/(m k)


        #length of the tank
        self.L_tank = 0.508 #m (~20in)
        #calculates tank volume
        self.V_tank = np.pi * (self.ro**2 - self.ri**2) * self.L_tank

        #percent of tank filled with liquid N2O
        self.fill_percent = 50 #% 

        #atmospheric temperature
        self.T_atm = 30 + 273.15 #K | atmospheric temperature
        self.P_atm = 101325 #Pa | atmospheric pressure 
        self.P_out = 6894800 #Pa | Pressure Inlet (Should eventually vary with source pressure)

        #timestep
        self.delta_t = 0 #s
        self.dV_vap_dt = 0 
        #Temperature of the Vapor part of wall
        self.T_wall_vap = self.T_atm #K | wall temperature at vapor region
        self.T_wall_liq = self.T_atm #K | wall temperature at liquid region
        self.T_lv_surf = 22 + 273.15 #K | inside temeprature at liquid vapor boundry
        self.T_liq = 20 + 273.15 #K | liquid temperature
        self.T_vap = 22 + 273.15 #K | vapor temperature
        
        
        N2O_mass = 15.0 #kg | Total mass of N2O
        fill_percent = 0 #% | percentage initially filled with liquid N2O
        t_stop = 600 #s () | total process time

        self.fill = True #TRUE: Tank is filling with N2O | FALSE: Tank is emptying N2O

        #eventually will chage m_dot to more accurately model pass flow rate
        if self.fill == True: # mass flow rate is negative with respect to the outlet (going in)
            self.M_dot = - N2O_mass / t_stop #kg/s |mass flow rate
            #All initial values will be 0 to reflect a perfectly empty tank (will add gaseous content later...)
            self.V_liq = 0
            self.V_vap = 0
            self.L_liq = 0 
            self.L_vap = 0
            self.A_liq_in = 0
            self.A_vap_in = 0
            self.A_liq_out = 0
            self.A_vap_out = 0
            self.m_wall_liq = 0
            self.m_wall_vap = 0
            self.V_wall_vap = 0
            self.m_liq = 0
            self.rho_liq = 0
            U_vap = 0
            U_liq = 0
            T_vap = 0
            T_liq = 0
        
        elif self.fill == False: #mass flow rate is positive with respect to the outlet (going out)
            self.M_dot = N2O_mass / t_stop #kg/s | mass flow rate
            #calculate initial length, area, and volume dimensions
            init_dim = self.Dim_Calcs((N2O_mass * fill_percent),PropsSI('D','T',self.T_liq,'Q',0,'N2O'))        
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
            self.V_wall_vap = np.pi * (self.ro**2 - self.ri**2) * self.L_vap

            self.m_liq = N2O_mass * fill_percent
            self.m_vap = (1 - fill_percent) * N2O_mass

            self.rho_liq = PropsSI('D','T',self.T_liq,'Q',0,'N2O')
            U_vap = self.m_vap * PropsSI('U','P' , self.P_out, 'T',self.T_vap, 'N2O')
            U_liq = self.m_liq * PropsSI('U','P' , self.P_out, 'T',self.T_liq, 'N2O')
            T_vap = self.T_atm - 2 #Change initial temp
            T_liq = self.T_atm - 3 #change initial temp

        ##Assuming all heat conduction equals 0 at start....
        self.z = [0, 0, 0, 0, 0, 0,\
             0, 0, 0, PropsSI('D','T',self.T_vap,'Q',1,'N2O'), self.L_vap, self.V_wall_vap ,\
             0, 0, 0, 0, 0, 0, self.T_atm-20, self.T_atm-20, self.m_liq, self.rho_liq, U_vap,\
                 U_liq, T_vap, T_liq ] 

        #self.x = odeint(self.ODE, self.z, self.tspan)
        self.t = 0 #s | time counter
        self.t_step = 0.1 #s |time step
        
        #######
        ## FOR TESTING PURPOSES ONLY
        #######
        
        self.t_plt = []
        self.plot_0 = []
        self.plot_1 = []
        self.plot_2 = []
        self.plot_3 = []
        self.plot_4 = []
        self.plot_5 = []
        self.plot_6 = []
        self.plot_7 = []
        self.plot_8 = []
        self.plot_9 = []
        self.plot_10 = []
        self.plot_11 = []
        self.plot_12 = []
        self.plot_13 = []
        self.plot_14 = []
        self.plot_15 = []
        self.plot_16 = []
        self.plot_17 = []
        self.plot_18 = []
        self.plot_19 = []
        self.plot_20 = []
        self.plot_21 = []

        self.Plot_store = []
        ##########################

        k = 1000

        for i in range(0,k):
            self.count_bar(i,k)

            ######
            # FOR TESTING PURPOSES ONLY
            ######
            self.t_plt.append(self.t)

            
            self.plot_0.append(self.z[0]) #Q_liq_surf
            self.plot_1.append(self.z[1]) #Q_surf_vap           
            self.plot_2.append(self.z[2]) #Q_wall_vap_out 
            self.plot_3.append(self.z[3]) #Q_wall_liq_out
            self.plot_4.append(self.z[4]) #m_evap
            self.plot_5.append(self.z[5]) #m_cond
            self.plot_6.append(self.z[6]) #m_vap
            self.plot_7.append(self.z[7]) #Q_in_vap
            self.plot_8.append(self.z[8]) #Q_in_liq
            self.plot_9.append(self.z[9]) #rho_vap 
            self.plot_10.append(self.z[10]) #L_vap                        
            self.plot_11.append(self.z[11]) #V_wall_vap
            self.plot_12.append(self.z[12]) #Q_wall_vap_in
            self.plot_13.append(self.z[13]) #Q_wall_liq_in 
            self.plot_14.append(self.z[14]) #Q_wall_vap_cond
            self.plot_15.append(self.z[15]) #Q_wall_liq_cond
            self.plot_16.append(self.z[16]) #m_wall_vap_in 
            self.plot_17.append(self.z[17]) #m_wall_liq_in
            self.plot_18.append(self.z[18]) #T_wall_vap
            self.plot_19.append(self.z[19]) #T_wall_liq
            self.plot_20.append(self.z[20]) #m_liq                       
            #self.plot_21.append(self.z[21]) #
            




            answer = self.ODE(self.z, self.t_step)
            #Initial conditions plus rate of change.... 
            self.z += np.dot(self.t_step, answer)
            self.t += self.t_step        

        self.Visualize()
        print(self.z)
        #print(new)

    def Visualize(self): # used to plot the data, toggle graphs on and off...
        #For Testing Purposes....


        plt.ylabel('Q_liq_surf')
        plt.plot(self.t_plt, self.plot_0) #Q_liq_surf
        plt.show()

        plt.ylabel('Q_surf_vap')    
        plt.plot(self.t_plt, self.plot_1, '-', lw=2) #Q_surf_vap           
        plt.show()

        plt.ylabel('Q_wall_vap_out')
        plt.plot(self.t_plt, self.plot_2, '-', lw=2) #Q_wall_vap_out 
        plt.show()

        plt.ylabel('Q_wall_liq_out')
        plt.plot(self.t_plt, self.plot_3) #Q_wall_liq_out
        plt.show()
        
        plt.ylabel('m_evap')
        plt.plot(self.t_plt, self.plot_4) #m_evap
        plt.show()        
        
        plt.ylabel('m_cond')
        plt.plot(self.t_plt, self.plot_5) #m_cond
        plt.show()
        
        plt.ylabel('m_vap')
        plt.plot(self.t_plt, self.plot_6) #m_vap
        plt.show()

        plt.ylabel('Q_in_vap')
        plt.plot(self.t_plt, self.plot_7) #Q_in_vap
        plt.show()

        plt.ylabel('Q_in_liq')
        plt.plot(self.t_plt, self.plot_8) #Q_in_liq
        plt.show()
        
        plt.ylabel('rho_vap')
        plt.plot(self.t_plt, self.plot_9) #rho_vap 
        plt.show()

        plt.ylabel('L_vap')
        plt.plot(self.t_plt, self.plot_10) #L_vap                        
        plt.show()
        
        plt.ylabel('V_wall_vap')
        plt.plot(self.t_plt, self.plot_11) #V_wall_vap
        plt.show()
        
        plt.ylabel('Q_wall_vap_in')
        plt.plot(self.t_plt, self.plot_12) #Q_wall_vap_in
        plt.show()
        
        plt.ylabel('Q_wall_liq_in')
        plt.plot(self.t_plt, self.plot_13) #Q_wall_liq_in 
        plt.show()
        
        plt.ylabel('Q_wall_vap_cond')
        plt.plot(self.t_plt, self.plot_14) #Q_wall_vap_cond
        plt.show()
        
        plt.ylabel('Q_wall_liq_cond')
        plt.plot(self.t_plt, self.plot_15) #Q_wall_liq_cond
        plt.show()
        
        plt.ylabel('m_wall_vap_in')
        plt.plot(self.t_plt, self.plot_16) #m_wall_vap_in 
        plt.show()
        
        plt.ylabel('m_wall_liq_in')
        plt.plot(self.t_plt, self.plot_17) #m_wall_liq_in
        plt.show()
        
        plt.ylabel('T_wall_vap')
        plt.plot(self.t_plt, self.plot_18) #T_wall_vap
        plt.show()
        
        plt.ylabel('T_wall_liq')
        plt.plot(self.t_plt, self.plot_19) #T_wall_liq
        plt.show()

        plt.ylabel('m_liq')
        plt.plot(self.t_plt, self.plot_20) #m_liq
        plt.show()


    def Dim_Calcs(self, m_liq, rho_liq):
        ###########
        ###NEED TO ADD CALCULATION FOR VOLUME OF WALL VAPOR AND LIQUID REGIONS
        ###########

        V_liq = m_liq / rho_liq #m^3 | volume of liquid in tank
        V_vap = self.V_tank - V_liq #m^3 | volume of vapor in tank
        if V_liq + V_vap > self.V_tank:
            print("Error in volume calcs...")
        
        L_liq = V_liq / (np.pi * self.ri**2) #m | Length of liquid region
        L_vap = V_vap / (np.pi * self.ri**2) #m | Length of vapor region
        if L_liq + L_vap > self.L_tank:
            print("Error in length calcs...")
        
        #Note::: Surface area does not include end caps; hollow tube assumed
        A_liq_in = 2 * np.pi * self.ri * L_liq #m^2 |inner surface area of tank liquid region
        A_vap_in = 2 * np.pi * self.ri * L_vap #m^2 |inner surface area of tank vapor region

        A_liq_out = 2 * np.pi * self.ro * L_liq #m^2 | outter surface area of tank liquid region
        A_vap_out = 2 * np.pi * self.ro * L_vap #m^2 | outter surface area of tank vapor region

        m_wall_liq = np.pi * (self.ro**2 - self.ri**2) * L_liq * self.rho_wall #kg | mass of wall in liquid region
        m_wall_vap = np.pi * (self.ro**2 - self.ri**2) * L_vap * self.rho_wall #kg | mass of wall in liquid region
        
        
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
        #Saturation Pressure
        P_sat = PropsSI('P','T',T,'Q',1,'N2O') #Pa
        #Molar Mass
        MM = PropsSI('M','T',T,'P',P,'N2O') #kg/mol
        #Compressibility Factor
        Z_vap = PropsSI('Z','T',T,'P',P,'N2O')
        #Gas Constant
        R_vap = PropsSI('GAS_CONSTANT','T',T,'P',P,'N2O') #J/mol/K 
        #The Following from... "Thermophysical properties of nitrous oxide" from IHS ESDU
        #Viscosity from EQ 4.9
        mu_vap = np.exp(3.3281 + (-1.18237) * (1/(T/309.59) - 1)**(1/3) + (-0.055155) * (1/(T/309.59) - 1)**(4/3)) #microNs/m^2
        mu_vap = mu_vap/ 1000000 #Ns/m^2
        #thermal conductivity
        k_vap =  np.exp(-7.0887 + (-0.276962) * (1-(T/309.59))**(-2/3) + 2.8872 * (1-(T/309.59))**(-1/3) + 16.6116 * \
                    (1-(T/309.59))**(1/3) + (-11.8221) * (1-(T/309.59))**(2/3)) #mW/m/K
        k_vap = k_vap/1000 #W/m/K
        return Cp_vap, rho_vap, mu_vap, k_vap, beta_vap, enthalpy_vap, enthalpy_vap_sat, P_sat, MM, Z_vap, R_vap


        #z = [self.T_wall_vap, 0, 0,   self.T_wall_liq, self.T_lv_surf, self.T_liq, self.T_vap, self.V_liq, self.V_vap,\
        #    self.L_liq, self.L_vap, self.A_liq_in, self.A_vap_in, self.A_liq_out, self.A_vap_out

        #z = [0, 0, 0, 0, 0     ]

    def Enthalpy_calc(self, P):
        ####
        # Should eventually vary with Pressure instead of density....
        ####
        h_evap = PropsSI('H', 'P', P, 'Q', 1, 'N2O') #J/kg | evaporation specific enthalpy
        h_cond = PropsSI('H', 'P', P, 'Q', 0, 'N2O') #J/kg | condensation specific enthalpy
        
        ######
        # h_out ...will be intresting to model with inward flow since it will be dependent on the source 
        # tank exit conditions and the losses which occur in the line. 
        ######
        
        h_out = PropsSI('H', 'P', P, 'Q', 0, 'N2O') #J/kg  | outlet specific enthalpy | Assumes liquid is going in/out....
        return h_evap, h_cond, h_out

    def ODE(self, z, delta_t):

        ###Need to rearrange initial conditions to match derivatives taken


        Q_liq_surf = z[0]        
        Q_surf_vap = z[1]        
        Q_wall_vap_out = z[2]  
        Q_wall_liq_out = z[3]        
        m_evap = z[4]
        m_cond = z[5]
        m_vap = z[6]
        Q_in_vap = z[7]
        Q_in_liq = z[8]
        rho_vap = z[9]
        L_vap = z[10]
        V_wall_vap = z[11]
        Q_wall_vap_in = z[12]
        Q_wall_liq_in = z[13]
        Q_wall_vap_cond = z[14]
        Q_wall_liq_cond = z[15]
        m_wall_vap_in = z[16]
        m_wall_liq_in = z[17]
        T_wall_vap = z[18]
        T_wall_liq = z[19]
        m_liq = z[20]
        rho_liq = z[21]
        U_vap = z[22]
        U_liq = z[23]
        T_vap = z[24]
        T_liq = z[25]
       

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

        dim = self.Dim_Calcs(m_liq, rho_liq)

        V_liq = dim[0]
        V_vap = dim[1]
        L_liq = dim[2]
        #L_vap = dim[3]
        A_liq_in = dim[4]
        A_vap_in = dim[5]
        A_liq_out = dim[6]
        A_vap_out = dim[7]
        m_wall_liq = dim[8]
        m_wall_vap = dim[9]

        P = 101325 * 800 
        self.dV_vap_dt = abs(V_vap) - abs(self.dV_vap_dt)


        #The film temperature between air and vapor wall region
        T_film_wall_vap = 0.5 * (z[0] + self.T_atm)
        #Finds air properties at film temp and atmospheric pressure
        Air = self.thermo_air(T_film_wall_vap, self.P_atm)
        Cp_air = Air[0]
        rho_air = Air[1]
        mu_air = Air[2]
        k_air = Air[3]
        beta_air = Air[4]


        T_film_liq_surf = 0.5 * (T_vap + T_liq)
        N2O_liq = self.Thermo_N2O_Liq(T_film_liq_surf, P)
        Cp_liq = N2O_liq[0]
        #rho_liq = N2O_liq[1]
        mu_liq = N2O_liq[2]
        k_liq = N2O_liq[3]
        beta_liq = N2O_liq[4]
        enthalpy_liq = N2O_liq[5]
        enthalpy_liq_sat = N2O_liq[6]

        N2O_vap = self.Thermo_N2O_Vap(T_film_liq_surf, P)
        Cp_vap = N2O_vap[0]
        #rho_vap = N2O_liq[1]
        mu_vap = N2O_vap[2]
        k_vap = N2O_vap[3]
        beta_vap = N2O_vap[4]
        #enthalpy_vap = N2O_vap[5]
        enthalpy_vap_sat = N2O_vap[6]
        P_sat = N2O_vap[7] +100
        MM_vap = N2O_vap[8]
        Z_vap = N2O_vap[9]
        R_vap = N2O_vap[10]


        N2O_latent_heat = enthalpy_vap_sat - enthalpy_liq_sat
        if L_liq == 0: #a hacky fix to avoid dividing by 0 in a filling model
            L_liq = 0.00000000001 #m |very hacky approach.......
        
        #dQ_liq_surf/dt
        EQ13 = (0.15 * ((Cp_liq * rho_liq**2 * g * beta_liq * abs(T_film_liq_surf - T_liq ) * L_liq**3 ) \
                / (mu_liq * k_liq ) )**(1/3) * (k_liq/L_liq)) * A_liq_in * (T_film_liq_surf - T_liq)
        
        #dQ_surf_vap/dt
        EQ14 = (0.15 * ((Cp_vap * rho_vap**2 * g * beta_vap * abs(T_film_liq_surf - T_vap) * L_vap**3 )   \
                / (mu_vap * k_vap ) )**(1/3) * (k_liq/L_vap)) * A_vap_in * (T_film_liq_surf - T_liq)
        
        #dQ_wall_vap_out/dt
        EQ2 = (0.021 * ((Cp_vap * rho_vap**2 * g * beta_vap * abs(T_wall_vap - T_vap) * L_vap**3  ) \
                / (mu_vap * k_vap ) )**(2/5) * (k_vap/L_vap)) * A_vap_in * (T_wall_vap - T_vap)

        #dQ_wall_liq_out/dt
        EQ7 = (0.021 * ((Cp_liq * rho_liq**2 * g * beta_liq * abs(T_wall_liq - T_liq) * L_liq**3  ) \
                / (mu_liq * k_liq ) )**(2/5) * (k_liq/L_liq)) * A_liq_in * (T_wall_liq - T_liq)

        #dm_evap/dt
        EQ12 = (EQ13 - EQ14) / (N2O_latent_heat + (enthalpy_liq_sat - enthalpy_liq) )
        

        #NEED TO CHECK CAN PRODUCE ERROR...
        
        #dm_cond/dt
        if P > P_sat:
            EQ15 = ( ( P - P_sat) * V_vap * MM_vap ) / (Z_vap * R_vap  * T_vap * delta_t )
        elif P <= P_sat: 
            EQ15 = 0

        #dm_vap
        EQ10 = EQ12 - EQ15

        #dQ_in_vap
        EQ24 = EQ2 + EQ14

        #dQ_in_liq
        EQ25 = EQ7 + EQ13


        #ERRORs.... Need a more accurate description of dV_vap_dt...
        
        #d_rho_vap/dt
        EQ18 = 1/V_vap * EQ10 - m_vap / (V_vap)**2 * self.dV_vap_dt
        
        #############
        ####NEED TO CHECK FOR UNCERTAINTY IN MEASURING THE HEIGHT.....
        #avoids having a 0 vapor density if the tank is completely empty, forces the Length of the vapor region to equal 0... 
        
        #dL_vap/dt
        if rho_vap == 0:
            EQ29 = 0
        else:
            EQ29 = 1/(np.pi * self.ri**2) * ((1/rho_vap) * EQ10 - (m_vap/rho_vap) * EQ18)

        ##### The rate of change of the vapor height needs to be checked....

        #dV_wall_vap/dt
        EQ31 = np.pi * (self.ro**2 - self.ri**2) * EQ29

        #dQ_wall_vap_in/dt
        EQ1 = (0.59 * (  (Cp_air * (rho_air )**2 * g * beta_air * abs(T_wall_vap - self.T_atm) * L_vap**3)  \
                / (mu_air * k_air ) )**(0.25) * (k_air/L_vap)) * A_vap_out * (T_wall_vap - self.T_atm)

        #dQ_wall_liq_in/dt
        EQ6 = (0.59 * (  (Cp_air * (rho_air )**2 * g * beta_air * abs(T_wall_liq - self.T_atm) * L_liq**3)  \
                / (mu_air * k_air ) )**(0.25) * (k_air/L_liq)) * A_liq_out * (T_wall_liq - self.T_atm)        


        #############
        # Note:
        # for EQ3 and EQ8 the sign of T_wall_liq - T_wall_vap may have to be reversed to that the heat transfer into a section
        # is proportional to te heat transfer out of a section maintaining an equvilancy.
        #############
        #dQ_wall_vap_cond/dt
        EQ3 = (self.k_wall * (T_wall_liq - T_wall_vap) * np.pi * (self.ro**2 - self.ri**2)) \
                / (0.5 * L_liq + 0.5 * L_vap)

        #dQ_wall_liq_cond/dt
        EQ8 = (self.k_wall * (T_wall_vap - T_wall_liq) * np.pi * (self.ro**2 - self.ri**2)) \
                / (0.5 * L_liq + 0.5 * L_vap)

        #dm_wall_vap_in/dt
        EQ4 = EQ31 * self.rho_wall
        
        #dm_wall_liq_in/dt 
        EQ9 = -EQ31 * self.rho_wall
        
        #dT_wall_vap/dt
        #DOUBLE CHECK T_wall_liq as T_w_in from EQs

        EQ0 = (EQ1 - EQ2 + EQ3 + EQ4 * self.C_wall * (T_wall_liq - T_wall_vap)) \
                / (m_wall_vap * self.C_wall)
        
        #dT_wall_vap/dt
        EQ5 = (EQ6 - EQ7 + EQ8 + EQ9 * self.C_wall * (T_wall_vap - T_wall_liq)) \
                / (m_wall_liq * self.C_wall)


        ###Note:
        # Change all M_dot later on...
        ###
        #dm_liq/dt
        EQ11 = -EQ12 + EQ15 - self.M_dot


        #####
        # Note:
        # Change V_liq to vary with time in a more elegant manner........
        #####
        #drho_liq/dt
        if V_liq == 0: #|temporarily fixes no liquid volume .....
            EQ21 = self.rho_liq
        else:
            EQ21 = 1/V_liq * EQ11 - (m_liq / V_liq**2) * (-self.dV_vap_dt)

        #####
        # assumes that the pressure at the exit is the same as the liquid pressure
        # which is the same as the vapor pressure
        #####

        H_vap_calc = self.Enthalpy_calc(self.P_out)

        #dU_vap/dt
        EQ17 = EQ12 * H_vap_calc[0] - EQ15 * H_vap_calc[1] - P * self.dV_vap_dt + EQ24
        
        #####
        # for EQ 20 check sign on "P *self.dV_vap_dt"... might be positive if using the vapor volume rate of change
        # 
        #####

        #dU_liq/dt
        EQ20 = self.M_dot * H_vap_calc[2] - EQ12 * H_vap_calc[0] + EQ15 * H_vap_calc[1] - P * self.dV_vap_dt + EQ25

        ########    
        # NEED TO ADD A TOTAL INTERNAL ENERGY COUNTER
        # FOR BOTH THE VAPOR AND LIQUID REGIONS TO THEN USE FOR THE 
        # SPECIFIC ENERGY CALCS....
        ########

        #dT_vap/dt
        EQ16 = 1 / Cp_vap * (1/m_vap * (EQ17 - EQ))

        #dT_liq/dit
        EQ19 = 
    
        #0 , 1, 2
        #h_evap, h_cond, h_out

        return EQ13, EQ14, EQ2, EQ7, EQ12, EQ15, EQ10, EQ24, EQ25, EQ18,\
                 EQ29, EQ31, EQ1, EQ6, EQ3, EQ8, EQ4, EQ9, EQ0, EQ5, EQ11, EQ21, \
                     EQ17, EQ20, EQ16, EQ19



    def count_bar(self,i,k):
        if i == 0:
            self.clear
            print("__________ | 0%")
        elif i/k == 0.1:
            self.clear
            print("#_________ | 10%")
        elif i/k == 0.2:
            self.clear
            print("##________ | 20%")        
        elif i/k == 0.3:
            self.clear
            print("###_______ | 30%")
        elif i/k == 0.4:
            self.clear
            print("####______ | 40%") 
        elif i/k == 0.5:
            self.clear
            print("#####_____ | 50%") 
        elif i/k == 0.6:
            self.clear
            print("######____ | 60%")
        elif i/k == 0.7:
            self.clear
            print("#######___ | 70%")        
        elif i/k == 0.8:
            self.clear
            print("########__ | 80%")
        elif i/k == 0.9:
            self.clear
            print("#########_ | 90%") 
        elif i/k >= 0.99999:
            self.clear
            print("########## | 100%")         

Model()