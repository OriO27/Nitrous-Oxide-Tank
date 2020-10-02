"""
OOrtiz
V: 1.3

Modeling Equations from 
"Review and Evalution of Models for 
Self-Pressurizing Propellant Tank Dynamics"
Zimmerman, J.E., Waxman, B.S., Cantwell, B.J.,
Zilliac, G.G. 49th Joint Propulsion Conference, 2013  

TO DO:
    Allow multiple material's to be selected as initial tank 
        properties
    Make tank outlet pressure (P_out) [or inlet from the filling
        perspective...] vary with source or destination pressure
        -Model with Bernoulli's EQ
    Model M_dot more accurately at tank outlet [inlet for filling]
    Fix Control Volume for liquid and vapor region's convective 
        surface area
    Add Air or an inert gas {N2} to the initially evacuated tank
    
NOTES:
    The print statements in the ODE method will be removed after testing
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
import os 

class Model(object):
    def __init__(self, Fill):
        self.clear = os.system('cls')
        #graviational acceleration
        global g
        g = 9.81 #m/s^2
        self.P_out = 6894800 #Pa | Pressure Inlet (Should eventually vary with source pressure)
        self.T_out = 298.15 #K | Temporarily bypasses the errors encountered during the filling process....
        self.Initial_Conditions(Fill)

        if Fill == True:
            self.z = self.Fill_Tank_IC()
        elif Fill == False:
            self.z = self.Empty_Tank_IC()

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
        self.plot_22 = []
        self.plot_23 = []
        self.plot_24 = []
        self.plot_25 = []
        self.plot_26 = []

        k = 2
        

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
            self.plot_21.append(self.z[21]) #rho_liq
            self.plot_22.append(self.z[22]) #U_vap
            self.plot_23.append(self.z[23]) #U_liq
            self.plot_24.append(self.z[24]) #T_vap
            self.plot_25.append(self.z[25]) #T_liq
            self.plot_25.append(self.z[26]) #T_liq



            answer = self.ODE(self.z, self.t_step)
            #Initial conditions plus rate of change.... 
            self.z += np.dot(self.t_step, answer)
            self.t += self.t_step        

        #self.Visualize()
        print(self.z)

    def Initial_Conditions(self, Fill):
        #Tank Outter (ro) and Inner (ri) radii || Diameter/2 
        self.ro = 0.1778/2 #m (7in)     
        self.ri = 0.1524/2 #m (6in)

        #Tank Properties (Material is Al for now)
        self.rho_wall = 2710 #kg/m3 | Density
        self.C_wall = 897 #J/ kg k | Heat Capacity
        self.k_wall = 150 #W/(m k) | Thermal Conductivity
        self.L_tank = 1.016 #m (40in) | Tank Length
        
        #Calculates Volume of Tank Material 
        self.V_tank_wall = np.pi * (self.ro**2 - self.ri**2) * self.L_tank #m^3

        #Calculates Tank Holding Volume (Inner Volume)
        self.V_tank_hold = np.pi * self.ri**2 * self.L_tank #m^3

        #Atmospheric Conditions 
        self.T_atm = 33 + 273.15 #K | Temperature
        self.P_atm = 101325 #Pa | Pressure 

        #Initial Temperature Conditions
        self.T_wall_vap = self.T_atm #K | wall temperature at vapor region
        self.T_wall_liq = self.T_atm #K | wall temperature at liquid region
        self.T_lv_surf = 22 + 273.15 #K | inside temeprature at liquid vapor boundry
        self.T_liq = 22 + 273.15 #K | liquid temperature
        self.T_vap = 22 + 273.15 #K | vapor temperature
        
        self.N2O_mass = 15.0 #kg | Total mass of N2O
        self.fill_percent = 50 #% | percentage initially filled with liquid N2O
        self.t_total = 600 #s () | Total process time...
        
        self.dV_vap_dt = 0 

    def Fill_Tank_IC(self):
        self.P_tank = 90000 #Pa | Initially Evacuated, will add Air
        self.M_dot = - self.N2O_mass / self.t_total #Mass Flow Rate at Entrance
        Tank_IC = [
            0, #Q_liq_surf     
            0, #Q_surf_vap     
            0, #Q_wall_vap_out
            0, #Q_wall_liq_out     
            0, #m_evap
            0, #m_cond
            0, #m_vap
            0, #Q_in_vap
            0, #Q_in_liq
            0, #rho_vap
            0, #L_vap
            0, #V_wall_vap
            0, #Q_wall_vap_in
            0, #Q_wall_liq_in
            0, #Q_wall_vap_cond
            0, #Q_wall_liq_cond
            0, #m_wall_vap_in
            0, #m_wall_liq_in
            0, #T_wall_vap
            0, #T_wall_liq
            0, #m_liq
            0, #rho_liq
            0, #U_vap
            0, #U_liq
            self.T_atm, #T_vap
            self.T_atm, #T_liq
            self.P_tank #P_tank
        ]
       
        return Tank_IC

    def Empty_Tank_IC(self):
        
        self.P_tank = 6850000 #Pa | Initial Tank Pressure
        self.M_dot = self.N2O_mass / self.t_total #kg/s | Mass Flow Rate at Exit
        
        H = PropsSI('H', 'T', self.T_liq, 'P', self.P_tank, 'N2O')
        H_liq = PropsSI('H', 'Q', 0, 'P', self.P_tank, 'N2O')

        if H <= H_liq: #Compressed Fluid Region
            m_liq = self.N2O_mass
            m_vap = 0
            rho_liq = self.N2O_mass/self.V_tank_hold
            rho_vap = 0
            U_liq = PropsSI('U','D', rho_liq, 'P', self.P_tank, 'N2O') #kg/m^3 |Liquid Density
            U_vap = 0

        Init_Dim = self.Dim_Calcs(m_liq, rho_liq)
        
        Tank_IC = [
            0, #Q_liq_surf     
            0, #Q_surf_vap     
            0, #Q_wall_vap_out
            0, #Q_wall_liq_out     
            0, #m_evap
            0, #m_cond
            m_vap, #m_vap
            0, #Q_in_vap
            0, #Q_in_liq
            rho_vap, #rho_vap
            Init_Dim[3], #L_vap
            Init_Dim[10], #V_wall_vap
            0, #Q_wall_vap_in
            0, #Q_wall_liq_in
            0, #Q_wall_vap_cond
            0, #Q_wall_liq_cond
            Init_Dim[9], #m_wall_vap_in
            Init_Dim[8], #m_wall_liq_in
            self.T_wall_vap, #T_wall_vap
            self.T_liq, #T_wall_liq
            m_liq, #m_liq
            rho_liq, #rho_liq
            U_vap, #U_vap
            U_liq, #U_liq
            self.T_vap, #T_vap
            self.T_liq, #T_liq
            self.P_tank #P_tank
        ]
        return Tank_IC

    def Visualize(self): # used to plot the data, toggle graphs on and off...
        #For Testing Purposes....


        plt.ylabel('Q_liq_surf')
        plt.xlabel('Time (s)')        
        plt.plot(self.t_plt, self.plot_0) #Q_liq_surf
        plt.show()

        plt.ylabel('Q_surf_vap')
        plt.xlabel('Time (s)')    
        plt.plot(self.t_plt, self.plot_1, '-', lw=2) #Q_surf_vap           
        plt.show()

        plt.ylabel('Q_wall_vap_out')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_2, '-', lw=2) #Q_wall_vap_out 
        plt.show()

        plt.ylabel('Q_wall_liq_out')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_3) #Q_wall_liq_out
        plt.show()
        
        plt.ylabel('m_evap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_4) #m_evap
        plt.show()        
        
        plt.ylabel('m_cond')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_5) #m_cond
        plt.show()
        
        plt.ylabel('m_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_6) #m_vap
        plt.show()

        plt.ylabel('Q_in_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_7) #Q_in_vap
        plt.show()

        plt.ylabel('Q_in_liq')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_8) #Q_in_liq
        plt.show()
        
        plt.ylabel('rho_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_9) #rho_vap 
        plt.show()

        plt.ylabel('L_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_10) #L_vap                        
        plt.show()
        
        plt.ylabel('V_wall_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_11) #V_wall_vap
        plt.show()
        
        plt.ylabel('Q_wall_vap_in')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_12) #Q_wall_vap_in
        plt.show()
        
        plt.ylabel('Q_wall_liq_in')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_13) #Q_wall_liq_in 
        plt.show()
        
        plt.ylabel('Q_wall_vap_cond')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_14) #Q_wall_vap_cond
        plt.show()
        
        plt.ylabel('Q_wall_liq_cond')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_15) #Q_wall_liq_cond
        plt.show()
        
        plt.ylabel('m_wall_vap_in')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_16) #m_wall_vap_in 
        plt.show()
        
        plt.ylabel('m_wall_liq_in')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_17) #m_wall_liq_in
        plt.show()
        
        plt.ylabel('T_wall_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_18) #T_wall_vap
        plt.show()
        
        plt.ylabel('T_wall_liq')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_19) #T_wall_liq
        plt.show()

        plt.ylabel('m_liq')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_20) #m_liq
        plt.show()

        plt.ylabel('rho_liq')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_21) #rho_liq
        plt.show()

        plt.ylabel('U_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_22) #U_vap
        plt.show()

        plt.ylabel('U_liq')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_23) #U_liq
        plt.show()

        plt.ylabel('T_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_24) #T_vap
        plt.show()

        plt.ylabel('T_vap')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_25) #T_liq
        plt.show()

        plt.ylabel('P')
        plt.xlabel('Time (s)')
        plt.plot(self.t_plt, self.plot_26) #P
        plt.show()

    def Dim_Calcs(self, m_liq, rho_liq):
        if rho_liq == 0:
            return 0,0,0,0,0,0,0,0,0,0,0,0

        V_liq = m_liq / rho_liq #m^3 | volume of liquid in tank
        V_vap = self.V_tank_hold - V_liq #m^3 | volume of vapor in tank
        
        if V_liq + V_vap > self.V_tank_hold:
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
            #6,7
        m_wall_liq = np.pi * (self.ro**2 - self.ri**2) * L_liq * self.rho_wall #kg | mass of wall in liquid region
        m_wall_vap = np.pi * (self.ro**2 - self.ri**2) * L_vap * self.rho_wall #kg | mass of wall in liquid region
            #8, 9
        V_wall_vap = np.pi * (self.ro**2 - self.ri**2) * L_vap #m^3 | volume of vapor region wall
        V_wall_liq = np.pi * (self.ro**2 - self.ri**2) * L_liq #m^3 | volume of liquid region wall
                #10, 11
        return V_liq, V_vap, L_liq, L_vap, A_liq_in, A_vap_in, A_liq_out,\
             A_vap_out, m_wall_liq, m_wall_vap, V_wall_vap, V_wall_liq


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

    def PDE(self, T, rho, m): #Used to approximate the partial derivative [del(u)/del(rho)]_T
        #+- 1% of the density
        if rho <= 0:
            return 0 #no approximate solution if density = 0

        rho_plus_1 = 1.01 * rho 
        rho_minus_1 = 0.99 * rho
        # using the high and low density values to calculate high and low 
        # internal energy values
        U_high = PropsSI('U', 'T', T, 'D', rho_plus_1, 'N2O')
        U_low = PropsSI('U', 'T', T, 'D', rho_minus_1, 'N2O')

        del_u = (U_high - U_low)/m #returns the specific change in internal energy
        del_rho = (rho_plus_1 - rho_minus_1)
        
        # approximate solution to the PDE
        del_u_rho = (del_u/del_rho) 
        
        return del_u_rho

    def ODE(self, z, delta_t):

        ###Need to rearrange initial conditions to match derivatives taken
        
        #############################################
        ## Unpack initial conditions to make the   ##
        ## equations easier to follow when working ##
        #############################################

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
        P = z[26]

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

        self.dV_vap_dt = abs(V_vap) - abs(self.dV_vap_dt)


        #The film temperature between air and vapor wall region
        T_film_wall_vap = 0.5 * (z[0] + self.T_atm)
        if T_film_wall_vap <= 183: #hacky approach to fix initial temperature starting at 0K due to evacuated tank... 
            T_film_wall_vap = self.T_out
        
        #Finds air properties at film temp and atmospheric pressure
        Air = self.thermo_air(T_film_wall_vap, self.P_atm)
        Cp_air = Air[0]
        rho_air = Air[1]
        mu_air = Air[2]
        k_air = Air[3]
        beta_air = Air[4]
        
        #TEST FOR AIR THERMO PROPERTIES 
        print('##############################')
        print('Cp_air = {}'.format(Air[0]))
        print('rho_air = {}'.format(Air[1]))
        print('mu_air = {}'.format(Air[2]))
        print('k_air = {}'.format(Air[3]))
        print('beta_air = {}'.format(Air[4]))





        T_film_liq_surf = 0.5 * (T_vap + T_liq)
        
        if T_film_liq_surf <= 183:
            T_film_liq_surf = self.T_out

        N2O_liq = self.Thermo_N2O_Liq(T_film_liq_surf, P)
        Cp_liq = N2O_liq[0]
        #rho_liq = N2O_liq[1]
        mu_liq = N2O_liq[2]
        k_liq = N2O_liq[3]
        beta_liq = N2O_liq[4]
        enthalpy_liq = N2O_liq[5]
        enthalpy_liq_sat = N2O_liq[6]
        print('##############################')
        print('Cp_liq = {}'.format(N2O_liq[0]))
        print('rho_liq = {}'.format(N2O_liq[1]))
        print('mu_liq = {}'.format(N2O_liq[2]))
        print('k_liq = {}'.format(N2O_liq[3]))
        print('beta_liq = {}'.format(N2O_liq[4]))
        print('enthalpy_liq = {}'.format(N2O_liq[5]))
        print('enthalpy_liq = {}'.format(N2O_liq[6]))
        print('##############################')
        N2O_vap = self.Thermo_N2O_Vap(T_film_liq_surf, P)
        Cp_vap = N2O_vap[0]
        #rho_vap = N2O_liq[1]
        mu_vap = N2O_vap[2]
        k_vap = N2O_vap[3]
        beta_vap = N2O_vap[4]
        #enthalpy_vap = N2O_vap[5]
        enthalpy_vap_sat = N2O_vap[6]
        P_sat = N2O_vap[7]
        MM_vap = N2O_vap[8]
        Z_vap = N2O_vap[9]
        R_vap = N2O_vap[10]

        print('Cp_vap = {}'.format(N2O_vap[0]))
        print('rho_vap = {}'.format(N2O_vap[1]))
        print('mu_vap = {}'.format(N2O_vap[2]))
        print('k_vap = {}'.format(N2O_vap[3]))
        print('beta_vap = {}'.format(N2O_vap[4]))
        print('enthalpy_vap = {}'.format(N2O_vap[5]))
        print('enthalpy_vap = {}'.format(N2O_vap[6]))
        print('P_sat = {}'.format(N2O_vap[7]))
        print('MM_vap = {}'.format(N2O_vap[8]))
        print('Z_vap = {}'.format(N2O_vap[9]))
        print('R_vap = {}'.format(N2O_vap[10]))





        #########################
        ## Unpacking compelete ##
        #########################

        N2O_latent_heat = enthalpy_vap_sat - enthalpy_liq_sat


        if L_liq == 0: #a hacky fix to avoid dividing by 0 in a filling model
            L_liq = 0.000000000001 #m |very hacky approach.......
        if L_vap == 0:
            L_vap = 0.000000000001 
        
        #######################################
        ## Note: All expressions are rounded ##
        ##   to 10 decimal places to avoid   ##
        ## potential RuntimeWarnings due to  ##
        ##  dividing by a very small number  ##
        #######################################

        #####################
        #   dQ_liq_surf/dt  #
        #####################
        EQ13 = (0.15 * ((Cp_liq * rho_liq**2 * g * beta_liq * abs(T_film_liq_surf - T_liq ) * L_liq**3 ) \
                / (mu_liq * k_liq ) )**(1/3) * (k_liq/L_liq)) * A_liq_in * (T_film_liq_surf - T_liq)
        EQ13 = round(EQ13, 10)
        print('dQ_liq_surf/dt = {}'.format(EQ13))

        #####################
        #   dQ_surf_vap/dt  #
        #####################
        EQ14 = (0.15 * ((Cp_vap * rho_vap**2 * g * beta_vap * abs(T_film_liq_surf - T_vap) * L_vap**3 )   \
                / (mu_vap * k_vap ) )**(1/3) * (k_liq/L_vap)) * A_vap_in * (T_film_liq_surf - T_liq)
        EQ14 = round(EQ14,10)
        print('dQ_surf_vap/dt = {}'.format(EQ14))

        ######################
        # dQ_wall_vap_out/dt #
        ######################
        EQ2 = (0.021 * ((Cp_vap * rho_vap**2 * g * beta_vap * abs(T_wall_vap - T_vap) * L_vap**3  ) \
                / (mu_vap * k_vap ) )**(2/5) * (k_vap/L_vap)) * A_vap_in * (T_wall_vap - T_vap)
        EQ2 = round(EQ2, 10)
        print('dQ_wall_vap_out/dt = {}'.format(EQ2))

        #####################
        #dQ_wall_liq_out/dt #
        #####################
        EQ7 = (0.021 * ((Cp_liq * rho_liq**2 * g * beta_liq * abs(T_wall_liq - T_liq) * L_liq**3  ) \
                / (mu_liq * k_liq ) )**(2/5) * (k_liq/L_liq)) * A_liq_in * (T_wall_liq - T_liq)
        EQ7 = round(EQ7, 10)
        print('dQ_wall_liq_out/dt = {}'.format(EQ7))

        ####################
        ##   dm_evap/dt   ##
        ####################
        EQ12 = (EQ13 - EQ14) / (N2O_latent_heat + (enthalpy_liq_sat - enthalpy_liq) )
        EQ12 = round(EQ12, 10)
        print('dm_evap/dt = {}'.format(EQ12))

        #NEED TO CHECK CAN PRODUCE ERROR...
        ####################
        #    dm_cond/dt    #
        ####################
        if P > P_sat:
            if T_vap <= 183:#attempts to fix an initial condition of 0K inside the tank....
                T_vap = self.T_out
            EQ15 = ( ( P - P_sat) * V_vap * MM_vap ) / (Z_vap * R_vap  * T_vap * delta_t )
            EQ15 = round(EQ15, 10)
        elif P <= P_sat: 
            EQ15 = 0

        print('dm_cond/dt = {}'.format(EQ15))

        #####################
        #       dm_vap      #
        #####################
        EQ10 = EQ12 - EQ15
        EQ10 = round(EQ10, 10)
        print('dm_vap = {}'.format(EQ10))

        #####################
        #      dQ_in_vap    #
        #####################
        EQ24 = EQ2 + EQ14
        EQ24 = round(EQ24, 10)
        print('dQ_in_vap = {}'.format(EQ24))

        ####################
        #     dQ_in_liq    #
        ####################
        EQ25 = EQ7 + EQ13
        EQ25 = round(EQ25, 10)
        print('dQ_in_liq = {}'.format(EQ25))

        #ERRORs.... Need a more accurate description of dV_vap_dt...
        #####################
        #     d_rho_vap/dt  #
        #####################
        if V_vap == 0:
            V_vap = 0.0000000000000001
        EQ18 = 1/V_vap * EQ10 - m_vap / (V_vap)**2 * self.dV_vap_dt
        EQ18 = round(EQ18, 10)
        print('d_rho_vap/dt = {}'.format(EQ18))

        #############
        ####NEED TO CHECK FOR UNCERTAINTY IN MEASURING THE HEIGHT.....
        #avoids having a 0 vapor density if the tank is completely empty, forces the Length of the vapor region to equal 0... 
        
        #####################
        #     dL_vap/dt     #
        #####################
        if rho_vap == 0:
            EQ29 = 0
        else:
            EQ29 = 1/(np.pi * self.ri**2) * ((1/rho_vap) * EQ10 - (m_vap/rho_vap) * EQ18)
            EQ29 = round(EQ29, 10)
        print('dL_vap/dt = {}'.format(EQ29))

        ##### The rate of change of the vapor height needs to be checked....

        #####################
        #   dV_wall_vap/dt  #
        #####################
        EQ31 = np.pi * (self.ro**2 - self.ri**2) * EQ29
        EQ31 = round(EQ31, 10)

        print('dV_wall_vap/dt = {}'.format(EQ31))

        #####################
        # dQ_wall_vap_in/dt #
        #####################
        EQ1 = (0.59 * (  (Cp_air * (rho_air )**2 * g * beta_air * abs(T_wall_vap - self.T_atm) * L_vap**3)  \
                / (mu_air * k_air ) )**(0.25) * (k_air/L_vap)) * A_vap_out * (T_wall_vap - self.T_atm)
        EQ1 = round(EQ1, 10)
        print('dQ_wall_vap_in/dt = {}'.format(EQ1))

        #####################
        # dQ_wall_liq_in/dt #
        #####################
        EQ6 = (0.59 * (  (Cp_air * (rho_air )**2 * g * beta_air * abs(T_wall_liq - self.T_atm) * L_liq**3)  \
                / (mu_air * k_air ) )**(0.25) * (k_air/L_liq)) * A_liq_out * (T_wall_liq - self.T_atm)        
        EQ6 = round(EQ6, 10)
        print('dQ_wall_liq_in/dt = {}'.format(EQ6))

        #############
        # Note:
        # for EQ3 and EQ8 the sign of T_wall_liq - T_wall_vap may have to be reversed to that the heat transfer into a section
        # is proportional to te heat transfer out of a section maintaining an equvilancy.
        #############

        #######################
        # dQ_wall_vap_cond/dt #
        #######################
        EQ3 = (self.k_wall * (T_wall_liq - T_wall_vap) * np.pi * (self.ro**2 - self.ri**2)) \
                / (0.5 * L_liq + 0.5 * L_vap)
        EQ3 = round(EQ3, 10)
        print('dQ_wall_vap_cond/dt = {}'.format(EQ3))

        #######################
        # dQ_wall_liq_cond/dt #
        #######################
        EQ8 = (self.k_wall * (T_wall_vap - T_wall_liq) * np.pi * (self.ro**2 - self.ri**2)) \
                / (0.5 * L_liq + 0.5 * L_vap)
        EQ8 = round(EQ8, 10)
        print('dQ_wall_liq_cond/dt = {}'.format(EQ8))

        #####################
        # dm_wall_vap_in/dt #
        #####################
        EQ4 = EQ31 * self.rho_wall
        EQ4 = round(EQ4, 10)
        print('dm_wall_vap_in/dt = {}'.format(EQ4))
        
        #####################
        # dm_wall_liq_in/dt #
        #####################
        EQ9 = -EQ31 * self.rho_wall
        EQ9 = round(EQ9, 10)
        print('dm_wall_liq_in/dt = {}'.format(EQ9))

        #DOUBLE CHECK T_wall_liq as T_w_in from EQs
        #####################
        #   dT_wall_vap/dt  #
        #####################
        try:
            EQ0 = (EQ1 - EQ2 + EQ3 + EQ4 * self.C_wall * (T_wall_liq - T_wall_vap)) \
                    / (m_wall_vap * self.C_wall)
        except:
            print('Division too small! Will approximate to 0!')
            EQ0 = 0
        else:
            EQ0 = round(EQ0, 10)
        print('dT_wall_vap/dt = {}'.format(EQ0))

        #####################
        #   dT_wall_vap/dt  #
        #####################
        try:
            EQ5 = (EQ6 - EQ7 + EQ8 + EQ9 * self.C_wall * (T_wall_vap - T_wall_liq)) \
                / (m_wall_liq * self.C_wall)
        except:
            print('Division too small! Will approximate to 0!')
            EQ5 = 0
        else:
            EQ5 = round(EQ5, 10)

        print('dT_wall_vap/dt = {}'.format(EQ5))

        ###Note:
        # Change all M_dot later on...
        ###
        #####################
        #     dm_liq/dt     #
        #####################
        EQ11 = -EQ12 + EQ15 - self.M_dot
        EQ11 = round(EQ11, 10)
        print('dm_liq/dt = {}'.format(EQ11))

        #####
        # Note:
        # Change V_liq to vary with time in a more elegant manner........
        #####

        #####################
        #    drho_liq/dt    #
        #####################
        if V_liq == 0: #|temporarily fixes no liquid volume .....
            EQ21 = rho_liq
        else:
            EQ21 = 1/V_liq * EQ11 - (m_liq / V_liq**2) * (-self.dV_vap_dt)
            EQ21 = round(EQ21, 10)

        print('drho_liq/dt = {}'.format(EQ21))

        #####
        # assumes that the pressure at the exit is the same as the liquid pressure
        # which is the same as the vapor pressure
        #####

        H_vap_calc = self.Enthalpy_calc(self.P_out)

        #####################
        #     dU_vap/dt     #
        #####################
        EQ17 = EQ12 * H_vap_calc[0] - EQ15 * H_vap_calc[1] - P * self.dV_vap_dt + EQ24
        EQ17 = round(EQ17,10)
        
        print('dU_vap/dt = {}'.format(EQ17))

        #####
        # for EQ 20 check sign on "P *self.dV_vap_dt"... might be positive if using the vapor volume rate of change
        # 
        #####

        #####################
        #     dU_liq/dt     #
        #####################
        EQ20 = self.M_dot * H_vap_calc[2] - EQ12 * H_vap_calc[0] + EQ15 * H_vap_calc[1] - P * self.dV_vap_dt + EQ25
        EQ20 = round(EQ20, 10)
        print('dU_liq/dt = {}'.format(EQ20))


        ########    
        # NEED TO ADD A TOTAL INTERNAL ENERGY COUNTER
        # FOR BOTH THE VAPOR AND LIQUID REGIONS TO THEN USE FOR THE 
        # SPECIFIC ENERGY CALCS....
        ########
        
        ####################
        #     dT_vap/dt    #
        ####################
        if m_vap == 0: # return 0 for temp if theres no N2O...
            EQ16 = 0
        else:
            EQ16 = (1 / Cp_vap) * (1/m_vap * (EQ17 - (U_vap / m_vap) * EQ10 ) - EQ18 * self.PDE(T_vap, rho_vap, m_vap) )
            EQ16 = round(EQ16, 10)
        print('dT_vap/dt = {}'.format(EQ16))

        ####################
        #     dT_liq/dt    #
        ####################
        if m_liq == 0: # return 0 for temp if theres no N2O...
            EQ19 = 0
        else:
            EQ19 = (1 / Cp_liq) * (1/m_liq * (EQ20 - (U_liq / m_liq) * EQ11 ) - EQ18 * self.PDE(T_liq, rho_liq, m_liq) )
            EQ19 = round(EQ19,10)
    
        print('dT_liq/dt = {}'.format(EQ19))

        return EQ13, EQ14, EQ2, EQ7, EQ12, EQ15, EQ10, EQ24, EQ25, EQ18,\
                 EQ29, EQ31, EQ1, EQ6, EQ3, EQ8, EQ4, EQ9, EQ0, EQ5, EQ11, EQ21, \
                     EQ17, EQ20, EQ16, EQ19, P



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

Model(True)