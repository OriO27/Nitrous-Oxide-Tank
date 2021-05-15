"""
OOrtiz
V: 3.7

Changes:
    N2O_Unsat treats compressed liquid viscosity and thermal conductivity (for T & P below the triple point) with the liquid equations from ESDU



TO DO:
Need to fix the pressure gradient at the pipe exit
Need to track total mass transfer through piping
Change draining initial conditons to be defined by T & P
Handle supercritical state from supply tank pressure (currently low to avoid supercritcal state)

"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp as e
from math import log as log
from scipy.optimize import root  # JS added
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI

class Main(object):
    def __init__(self):
        """
        t_step         :: Changes the time interval between integration solutions

        k              :: Is the total number of integration iterations
                                Note: the integration will automatically stop if the mass falls
                                below 0kg or above the initial tank mass
        N2O_total_mass :: Total N2O mass to be filled or emptied
        P_atm          :: Atmospheric pressure
        Fill           :: Set to True if modeling a filling process or False if modeling an emptying process
        end_loop       :: Should be False by default, used to terminate the iteration under certain circumstances
        z              :: Holds the data at each time step before appending to 'plot_x' lists
        plot_'x'       :: Lists that hold N2O state data

        Tank_Geometry():: Method where the N2O tank geometry (Length, Radius, Thickness) can be changed
        Pipe_info()    :: Method where the N2O pipe geometry (Total Length, Diameter) and total loss coefficients due to valves and turns, can be changed

        Tank_Initial() :: Method where the N2O and atmospheric temperatures can be changed
        Model()        :: Method that carries out the integration
        Visualize()    :: Method that plots the N2O state data
        """

        # Global Variables
        global g, t_step
        g = 9.81 #m/s^2
        t_step = 0.01 #s

        k = 30000
        self.N2O_total_mass = 15 #kg  
        self.t = 0.0
        self.P_atm = float(101325) #Pa | atmospheric pressure
        self.Fill = True #True for Fill | False for Empty

        self.end_loop = False
        self.Tank_Geometry()
        self.Pipe_Info()
        self.z = self.Tank_initial(self.Fill)

         # JS 3/11/2021 Added initial state
        self.t_plt = [0]     # Time
        self.plot_0 = [self.z[0]]    # N2O Temp
        self.plot_1 = [self.z[1]]    # U_N2O
        self.plot_2 = [self.z[2]]    # Wall Temp
        self.plot_3 = [self.z[3]]    # N2O Mass
        self.plot_4 = [self.z[4]]    # Tank Pressure
        self.plot_5 = [self.z[5]]    # Quality
        self.plot_6 = [self.z[6]]    # P_out

        try:
            self.Model(k)
            print("Total Process Time = {:.2f} s".format(self.t_plt[-1]))
            print("Final Tank Wall Temperature = {:.2f} K".format(self.plot_2[-1]))           
            print("Final N2O Temperature = {:.2f} K".format(self.plot_0[-1]))
            print("Final N2O Pressure = {:.3f} Pa".format(self.plot_4[-1]))           
            print("Final N2O Mass = {:.8f} kg".format(self.plot_3[-1]))  
            print("Final N2O Quality = {:.2f}".format(self.plot_5[-1]))              
            print("Final N2O Internal Energy = {:.2f} J".format(self.plot_1[-1]))
            self.Visualize()

        except:
            print('Error.... @ Line 79')
            print("Total Process Time = {:.2f} s".format(self.t_plt[-1]))
            print("Final Tank Wall Temperature = {:.2f} K".format(self.plot_2[-1]))           
            print("Final N2O Temperature = {:.2f} K".format(self.plot_0[-1]))
            print("Final N2O Pressure = {:.3f} Pa".format(self.plot_4[-1]))           
            print("Final N2O Mass = {:.8f} kg".format(self.plot_3[-1]))  
            print("Final N2O Quality = {:.2f}".format(self.plot_5[-1]))              
            print("Final N2O Internal Energy = {:.2f} J".format(self.plot_1[-1]))              
            self.Visualize()

    def Model(self, k):
        """
        This method loops k times in order to solve for N2O state data.
        Calls EQM_Tank_Solve() method, where the ODE's are solved.

        self.z is used as a temporary information holder.

        The info returned from EQM_Tank_Solve() method will be used as the intial
        condition of the Pipe_Solve() method to solve for the pipe state which will then
        influence the tank exit mass flow.

            Note:
            If self.end_loop is changed to True within one of the nested methods the solver will end
            If EQM_Tank_Solve() method returns 'next' the solver will skip appending any data and start
            on the next iteration
        """

        for solver in range(0,k):
            print("\n\n")
            print(solver) #counter check
            
            EQ_Solver = self.EQM_Tank_Solve(self.z)
            self.t += t_step
            if self.end_loop == True:
                print("Loop break @ solver start LINE 106")
                break
            elif EQ_Solver == 'next':
                continue

    ###################### JS major Change Summary #####################
    # PipeSolver calculates a huge pressure drop, such that the pressure is
    # negative at the exit of the tubing. So, you have an unphysically high
    # mass flow rate, that cannot exist with a tube of this length and diameter.
    # To deal with this, I think that your method of solving for the mass flow
    # shouldn't be used in filling mode. I'm trying a different method based on
    # The Bernoulli equation here.


            self.z = [
                    EQ_Solver[0],
                    EQ_Solver[1],
                    EQ_Solver[2],
                    EQ_Solver[3],
                    EQ_Solver[4],
                    EQ_Solver[5],
                    self.z[6]
                    ]

            self.t_plt.append(self.t)
            self.plot_0.append(self.z[0]) # N2O Temp
            self.plot_1.append(self.z[1]) # U_N2O
            self.plot_2.append(self.z[2]) # Wall Temp
            self.plot_3.append(self.z[3]) # N2O Mass
            self.plot_4.append(self.z[4]) # Tank Pressure
            self.plot_5.append(self.z[5]) # Quality
            self.plot_6.append(self.z[6]) # P_out...

            print('Write Try')
            Fill_Data = open("C:/Users/m_i_d/Desktop/Fill_Data.txt", "a+")
            print('write open')
            Fill_Data.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.t, self.z[0], self.z[1], self.z[2], self.z[3], self.z[4], self.z[5], self.z[6]))
            print('write down')
            Fill_Data.close()
            print('write close')
            # Ends the iteration if the mass falls below 0kg or it goes above 0.05% of the initial tank mass
            if self.z[3] <= 0 or self.z[3] >= 1.0005*self.N2O_total_mass:
                print("Loop Break\nMass out of range LINE 141")
                break
            elif self.t > 0.3:
                if self.plot_3[-1] - self.plot_3[-2] <= 0:
                    print("Change in Tank Mass Less Than 0.000006kg/min\nProgram Terminated")
                    break

    def Tank_initial(self, Fill):
        """
        Sets the initial temperature conditions for the tank and atmosphere

        Assumes atmospheric pressure is 101kPa

        Solves for the initial internal energy of the tank

        Sets the initial quality and tank pressures -- depends on selected filling or emptying mode
        """

        self.T_N2O = 294 #K | wall temperature at liquid region  # JS 3/11/2021 I made it colder to avoid supercritical T
        self.T_atm = 297 #K
        self.T_wall = 297

        if Fill == False: #Sets Initial Mass and Quality if tank is in Emptying Mode
            self.N2O_init_mass = self.N2O_total_mass #kg
            X = 0
            P_Tank = float(PropsSI('P', 'Q', X, 'T', self.T_N2O, 'N2O'))
            P_out = self.P_atm #Pa
            U_init = float(PropsSI('U', 'Q', X, 'T', self.T_N2O, 'N2O') * self.N2O_init_mass )  # JS 3/11/2021
        elif Fill == True: #Initial Mass and Quality if tank is in Filling Mode
            # self.N2O_init_mass = 0.001 #kg # JS 3/11/2021
            # X = 1 # JS 3/11/2021
            X = 2  # should be superheated, initially.
            P_Tank = self.P_atm #Pa | if Filling Tank Assumed to be at P_atm
            self.N2O_init_mass = PropsSI('D', 'T', self.T_N2O, 'P', P_Tank, 'N2O') * self.V_tank_inside # JS 3/11/2021
            # P_out = float(PropsSI('P', 'Q', X, 'T', self.T_N2O, 'N2O'))# JS 3/11/2021

            P_out = 68 * self.P_atm # This is the pressure of the supply tank, right? We can find this, I just made a guess  # JS 3/11/2021

            U_init = float(PropsSI('U', 'T', self.T_N2O, 'P', P_Tank, 'N2O') * self.N2O_init_mass)  # JS 3/11/2021

        # U_init = float(PropsSI('U', 'Q', X, 'T', self.T_N2O, 'N2O') * self.N2O_init_mass )  # JS 3/11/2021

        return self.T_N2O, U_init, self.T_wall, self.N2O_init_mass, P_Tank, X, P_out

    def Tank_Geometry(self):
        """
        Sets tank geometry, calls Tank_Material() method which enables density, heat capacity, and thermal conductivity of
            selected material. Made a separate method to eventaully allow for the slection of other materials.
        The tank exit diameter is used for the mass flow rate into/outof the tank and outof/into the pipe.
        """

        self.Tank_Material()
        self.L_tank = 1.524 # m (50in)   | Tank Length
        self.ri = 0.1524/2 # m (6in)    | Inner tank radius
        self.ro = 0.1778/2 # m (7in)    | Outter tank radius
        self.V_tank_inside = np.pi * self.ri**2 * self.L_tank # m^3 | Inner tank volume
        self.V_tank_wall = np.pi * (self.ro**2 - self.ri**2) * self.L_tank # m^3 | volume of the tank wall; used for heat conduction eqs
        self.A_tank_in = 2 * np.pi * self.ri * self.L_tank # m^2    | Inner surface area
        self.A_tank_out = 2 * np.pi * self.ro * self.L_tank # m^2   | Uutter surface area
        self.m_wall = np.pi * (self.ro**2 - self.ri**2) * self.L_tank * rho_wall # kg   | Tank wall mass
        self.R_cyl = float(np.log(self.ro/self.ri)/(2*np.pi*self.L_tank*k_wall))

        D_out = 0.00635 #m | Exit Diameter 1/4"
        self.A_tank_exit = np.pi * (D_out/2)**2 # m^2   | Tank exit area

    def Pipe_Info(self):
        """
        Method is used to set the pipe geometry/characteristics (length, diamter, area, roughness ) and loss coefficients

        k       ::  A list that holds individual loss coefficients

        K_pipe  ::  The sum of all the loss coefficients to be used in Pipe_Solve() method
        """
        self.D_pipe = 0.00635 #m | Pipe Diameter 1/4"
        self.A_pipe = np.pi/4 * self.D_pipe**2 #m^2 | Pipe Cross-Sectional Area
        self.rough = 0.000045 #m | Commercial Steel
        self.L_pipe = 9.144 # 3.048 #m | Total length of pipe 10ft
        k = [1, .5, 2.3, 2.3, 2.3, 70, 0.75, 2, 10, 6, 6, 6 ] #Loss coefficients due to bends, valves, etc. assume constant Velocity but inaccurate

        self.K_pipe = sum([losses for losses in k])

    def Tank_Material(self):
        """
        Method holds and returns tank material density, heat capacity, and thermal conductivity.
        """
        #Aluminum
        global rho_wall, C_wall, k_wall
        rho_wall = 2710.0 #kg/m3 | Density
        C_wall = 887.0 #J/( kg K) | Heat Capacity
        k_wall = 150.0 #W/(m K) | Thermal Conductivity

    def Air_prop(self, T):
        """
        When called, the method returns the heat capacity, density, viscosity, thermal conductivity, and expansion coefficients of air.
            Input temperature is the film temperature between the atmosphere and tank wall
        """
        Cp = float(PropsSI('C', 'T', T, 'P', self.P_atm, 'AIR'))   # J/kg/K    | Heat Capacity
        rho = float(PropsSI('D', 'T', T, 'P', self.P_atm, 'AIR'))  # kg/m^3    | Density
        k = float(PropsSI('L', 'T', T, 'P', self.P_atm, 'AIR'))    # 1/K       | Thermal Conductivity
        mu = float(PropsSI('V', 'T', T, 'P', self.P_atm, 'AIR'))   # Pa*s      | Viscosity
        beta = float(PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', T, 'P', self.P_atm, 'AIR')) # 1/K      | Expansion Coefficient

        return Cp, rho, mu, k, beta

    def N2O_prop(self, T, X):
        """
        Method calculates and returns thermodynamic properties of N2O at the input temperature and quality.
        The maximum allowed temperature input is the critical temperature, if any temperature

        Dynamic Viscosity and Thermal Conductivity From:
        "Thermophysical properties of nitrous oxide" IHS ESDU
        CoolProp does solve for these two properties for N2O

        Viscosity and thermal conductivity are solved when the quality is 0 and 1
            if the input quality is between 0 and 1 it is used to calculate the property
        """
        if not 0 <= X <= 1:
            raise ValueError('N2O_prop is only implemented for saturated mixtures')

        T_r = T/309.52 #reduced Temp

        Cp = float(PropsSI('C', 'T', T, 'Q', X, 'N2O'))  # J/kg/K      | Heat Capacity
        
        rho = float(PropsSI('D', 'T', T, 'Q', X, 'N2O')) # kg/m^3      | Density
        beta = float(PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', T, 'Q', X, 'N2O')) # 1/K   | Expansion Coefficient

        # Saturated Thermal Conductivity
        f1, f2, f3, f4 = 72.35, 1.5, -3.5, 4.5
        g1, g2, g3, g4, g5 = -7.08870, -0.276962, 2.88672, 16.6116, -11.8221
        # EQ 4.11 From ESDU | W/mK
        k_l = (f1 * (1 + f2 * (1-T_r)**(1/3) + f3*(1-T_r)**(2/3) + f4*(1-T_r)) )* 10**(-3)

        # EQ 4.12 From ESDU | W/mK
        lnk_v = g1 + g2*(1-T_r)**(-2/3) + g3*(1-T_r)**(-1/3) + g4*(1-T_r)**(1/3) + g5*(1-T_r)**(2/3)
        k_v = e(lnk_v)* 10**(-3)

        # Saturated Dynamic Viscosity
        h1, h2, h3, h4 = 1.6089, 2.0439, 5.24, 0.0293423
        i1, i2, i3 = 3.3281, -1.18237, -0.055155
        theta = (309.52 - h3) / (T - h3)

        # EQ 4.9 | kg/(m s)
        mu_l = (h4 * e(h1*(theta - 1)**(1/3) + h2*(theta - 1)**(4/3)) )* 10**(-3)

        # EQ 4.10 | kg/(m s)
        lnmu_v = i1 + i2*(1/T_r - 1)**(1/3) + i3*(1/T_r - 1)**(4/3)
        mu_v = e(lnmu_v) * 10**(-3)

        mu = (X * (mu_v - mu_l) + mu_l)
        k = (X * (k_v - k_l) + k_l)

        return Cp, rho, mu, k, beta

    def N2O_Unsat(self, T, P, m):
        """
        Method calculates and returns thermodynamic properties of N2O at the input temperature and pressure.
        Only valid for unsaturated conditions between the triple point and critical points.
        The maximum allowed temperature input is the critical temperature, if any temperature
        """
        try:
            phase_Check = CP.PhaseSI('P', P, 'T', T, 'N2O')
            print('Phase :: {}'.format(phase_Check))
        except:
            phase_Check = 'compressed_liquid'
            print('Phase :: {}'.format(phase_Check))       
            
        if 182.33 < T < 309.52:
            T_r = T/309.52 #reduced Temp
            # Thermal Conductivity
            f1, f2, f3, f4 = 72.35, 1.5, -3.5, 4.5
            g1, g2, g3, g4, g5 = -7.08870, -0.276962, 2.88672, 16.6116, -11.822
            # Dynamic Viscosity
            h1, h2, h3, h4 = 1.6089, 2.0439, 5.24, 0.0293423
            i1, i2, i3 = 3.3281, -1.18237, -0.055155
            theta = (309.52 - h3) / (T - h3)

            print('N2O_Unsat attempt')
            Cp = float(PropsSI('C', 'T', T, 'P', P, 'N2O'))  # J/kg/K      | Heat Capacity
            rho = float(PropsSI('D', 'T', T, 'P', P, 'N2O')) # kg/m^3      | Density
            beta = float(PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', T, 'P', P, 'N2O')) # 1/K   | Expansion Coefficient
            print("N2O_Unsat Cp | rho | beta found")
            if phase_Check == 'liquid':
                # EQ 4.11 From ESDU | W/mK
                k = (f1 * (1 + f2 * (1-T_r)**(1/3) + f3*(1-T_r)**(2/3) + f4*(1-T_r)) )* 10**(-3)
                # EQ 4.9 | kg/(m s)
                mu = (h4 * e(h1*(theta - 1)**(1/3) + h2*(theta - 1)**(4/3)) )* 10**(-3)
            elif phase_Check == 'gas':
                # EQ 4.12 From ESDU | W/mK
                lnk_v = g1 + g2*(1-T_r)**(-2/3) + g3*(1-T_r)**(-1/3) + g4*(1-T_r)**(1/3) + g5*(1-T_r)**(2/3)
                k = e(lnk_v)* 10**(-3)               
                # EQ 4.10 | kg/(m s)
                lnmu_v = i1 + i2*(1/T_r - 1)**(1/3) + i3*(1/T_r - 1)**(4/3)
                mu = e(lnmu_v) * 10**(-3)
        elif T >= 309.52:
            rho = m/self.V_tank_inside
            _, _, _, _, _, Cp, P = self.Helmholtz_alpha_eqs(T, (rho))
            beta = 1/T
            if 309.52 <= T <= 398.15:
                print("\nSupercritical Temperature = {}".format(T))
                print("Supercritical Pressure = {}".format(P))

                k = self.UnSat_Thermal_Cond(T, P)
                mu = self.Unsat_Visc(T, P)


        return Cp, rho, mu, k, beta

    def Mass_flow_JS(self, Ps, rho, mu):
        """Calculate the mass flow rate

        Governing Equation:
        (P2 - P1) = rho * g * head_loss

        Parameters
        ----------
        Ps : tuple of length 2
            Two floats indicating the pressure in the tank and at the far end
            of the tube. Units are Pascals
        rho : float
            Density of the liquid in the tube, kg/m^3
        mu : float
            Viscosity of the liquid in the tube, kg/(m s)

        Returns
        -------
        mdot : float
            Mass flow rate, kg/s
        """
        if Ps[0] > Ps[1]:   # JS 3/11/2021
            # There is a check valve, so flow cannot go backwards.
            return 0   # JS 3/11/2021
        h_loss = (Ps[1] - Ps[0]) / (rho * g)   # JS 3/11/2021

        def h_loss_func(V, rho, D, mu, L, K, rough, h_loss_expected):
            """Return the difference between head loss and expected head loss."""
            Re = (rho * V * D)/ mu

            if Re <= 2000:
                f = 64/Re
            else:
                # Haaland Equation
                f = ( 1 / \
                        (-1.8*log( ( (rough/D)/(3.7) )**1.11 + \
                            6.9/Re)))**2

            return V**2/(2 * g) * (K + f * L / D) - h_loss_expected

        # Solve the function above for its root
        args = (rho, self.D_pipe, mu, self.L_pipe, self.K_pipe, self.rough, h_loss)
        sol = root(h_loss_func, 1, args=args)
        velocity = sol.x[0]
        # if self.Fill:   # JS 3/11/2021
        #     sign = 1   # JS 3/11/2021
        # else:   # JS 3/11/2021
        #     sign = -1   # JS 3/11/2021

        return velocity * self.A_pipe * rho# * sign   # JS 3/11/2021

    def Mass_flux(self, A, P1, h1, rho1, P2, h2, rho2):
        """
        Method calcualtes the mass flux at the tank exit/entry based on the current tank pressure and the pipe pressure
            The mass flow direction is determined based on the pressure gradient, if else used to avoid possible
            negative square roots
        P1, h1, rho1    :: Pressure, enthalpy, & density of N2O in tank

        P2, h2, rho2    :: Pressure, enthalpy, & density of N2O at pipe entrance/exit

        Cd              :: Discharge coefficient of tank/pipe entrance/exit
            'end_loop' variable set to 'True' when in emptying mode and the tank pressure falls
            below the atmospheric pressure
        """

        Cd = 0.2

        if self.Fill == False:
            if P1 <= self.P_atm: #returns true to end loop
                self.end_loop = True
                print("P1 = {:.3e}".format(P1))
                print("P2 = {:.3e}".format(P2))
                print('Exit @ Mass Flow || P2 <= P_atm | L442')
                return self.end_loop

        # Simplified by JS
        G_SPI = Cd * (2 * rho1 * abs(P1-P2))**0.5  # Units are kg/(m^2 s)
        G_HEM = Cd * rho2 * (2 * abs(h1 - h2))**0.5  # Units are kg/(m^2 s)

        P1_sat = P1  # This isn't always true.......JS
        k = ((P1-P2)/(P1_sat - P2))**(.5)  # With the line above, this is always = 1.......JS
        G = float((k * G_SPI + G_HEM)/( 1 + k))

        if P1 >= P2: # Sets mass flow direction based on Pressure gradient
            m_dot = -G * A
        elif P1 < P2:
            m_dot = G * A

        #Used as a temporary check
        print("P1 = {:.1f}".format(P1))
        print("P2 = {:.1f}".format(P2))
        print('k = {}'.format(k))
        print('G_SPI = {:4e}'.format(G_SPI))
        print('G_HEM = {:4e}'.format(G_HEM))
        print('G = {:4e}'.format(G))
        print('m_dot = {:4e}'.format(m_dot))

        return m_dot

    def Helmholtz_alpha_eqs(self, T, rho):
        """
        Uses the Helmholtz fundamental equations to find the N2O state
        From: Short Fundamental Equations of State for 20 Industrial Fluids
        By: E.W. Lemmon & R. Span
        """
        T_c = 309.52
        rho_c = 452
        R_m  = 8.31446261815324 # J/(mol K) Molar Gas Constant
        m_w = (44.0128/1000) #kg/mol N2O Molecular Weight 
        R = R_m/m_w #J/(kg K) N2O Specific Gas Constant

        delta = rho/rho_c
        tau = T / T_c

        #Alpha_r constants
        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12 = \
            0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032, \
                -0.0036985, -0.23263, -0.00042859, -0.042810, -0.023038

        #Alpha_0 constants
        v1, u1, v2, u2, v3, u3 = 2.1769, 879.0, 1.6145, 2372.0, 0.48393, 5447.0
        c0, a1, a2 = 3.5, -4.4262736272, 4.3120475243

        alpha_r = n1 * delta * (tau**0.25) +\
                    n2 * delta * (tau**1.25) +\
                    n3 * delta * (tau**1.5) + \
                    n4 * (delta**3) * (tau**0.25) + \
                    n5 * (delta**7) * (tau**0.875) + \
                    n6 * delta * (tau**2.375) * (e(-delta)) + \
                    n7 * (delta**2) * (tau**2) * (e(-delta)) + \
                    n8 * (delta**5) * (tau**2.125) * (e(-delta)) + \
                    n9 * (delta) * (tau**3.5) * (e(-delta**2)) + \
                    n10 * delta * (tau**6.5) * (e(-delta**2)) + \
                    n11 * (delta**4) * (tau**4.75) * e(-delta**2) + \
                    n12 * (delta**2) * (tau**12.5) * e(-delta**3)

        # Partial alpha(r) with respect to partial delta 
        dalpha_r_ddelta = 7*delta**6*n5*tau**0.875 - 2*delta**5*n11*tau**4.75*e(-delta**2) -\
                    delta**5*n8*tau**2.125*e(-delta) - 3*delta**4*n12*tau**12.5*e(-delta**3) +\
                    5*delta**4*n8*tau**2.125*e(-delta) + 4*delta**3*n11*tau**4.75*e(-delta**2) -\
                    2*delta**2*n10*tau**6.5*e(-delta**2) + 3*delta**2*n4*tau**0.25 - delta**2*n7*tau**2*e(-delta) -\
                    2*delta**2*n9*tau**3.5*e(-delta**2) + 2*delta*n12*tau**12.5*e(-delta**3) - delta*n6*tau**2.375*e(-delta) +\
                    2*delta*n7*tau**2*e(-delta) + n1*tau**0.25 + n10*tau**6.5*e(-delta**2) + n2*tau**1.25 + n3*tau**1.5 +\
                    n6*tau**2.375*e(-delta) + n9*tau**3.5*e(-delta**2) 
            
        # 2nd Partial alpha(r) with respect to partial delta 
        d2alpha_r_ddelta2 = 4*delta**6*n11*tau**4.75*e(-delta**2) + 9*delta**6*n12*tau**12.5*e(-delta**3) +\
                    42*delta**5*n5*tau**0.875 + delta**5*n8*tau**2.125*e(-delta) - 18*delta**4*n11*tau**4.75*e(-delta**2) -\
                    10*delta**4*n8*tau**2.125*e(-delta) + 4*delta**3*n10*tau**6.5*e(-delta**2) - 18*delta**3*n12*tau**12.5*e(-delta**3) +\
                    20*delta**3*n8*tau**2.125*e(-delta) + 4*delta**3*n9*tau**3.5*e(-delta**2) + 12*delta**2*n11*tau**4.75*e(-delta**2) +\
                    delta**2*n7*tau**2*e(-delta) - 6*delta*n10*tau**6.5*e(-delta**2) + 6*delta*n4*tau**0.25 + delta*n6*tau**2.375*e(-delta) -\
                    4*delta*n7*tau**2*e(-delta) - 6*delta*n9*tau**3.5*e(-delta**2) + 2*n12*tau**12.5*e(-delta**3) - 2*n6*tau**2.375*e(-delta) + 2*n7*tau**2*e(-delta)

        # Partial alpha(r) with respect to partial tau
        dalpha_r_dtau = 0.875*delta**7*n5*tau**(-0.125) + 2.125*delta**5*n8*tau**1.125*e(-delta) +\
                    4.75*delta**4*n11*tau**3.75*e(-delta**2) + 0.25*delta**3*n4*tau**(-0.75) +\
                    12.5*delta**2*n12*tau**11.5*e(-delta**3) + 2*delta**2*n7*tau*e(-delta) +\
                    0.25*delta*n1*tau**(-0.75) + 6.5*delta*n10*tau**5.5*e(-delta**2) + 1.25*delta*n2*tau**0.25 +\
                    1.5*delta*n3*tau**0.5 + 2.375*delta*n6*tau**1.375*e(-delta) + 3.5*delta*n9*tau**2.5*e(-delta**2) 

        # 2nd Partial alpha(r) with respect to partial tau
        d2alpha_r_dtau2 = delta*(-0.109375*delta**6*n5*tau**(-1.125) + 2.390625*delta**4*n8*tau**0.125*e(-delta) +\
                    17.8125*delta**3*n11*tau**2.75*e(-delta**2) - 0.1875*delta**2*n4*tau**(-1.75) +\
                    143.75*delta*n12*tau**10.5*e(-delta**3) + 2*delta*n7*e(-delta) - 0.1875*n1*tau**(-1.75) +\
                    35.75*n10*tau**4.5*e(-delta**2) + 0.3125*n2*tau**(-0.75) + 0.75*n3*tau**(-0.5) +\
                    3.265625*n6*tau**0.375*e(-delta) + 8.75*n9*tau**1.5*e(-delta**2))

        d2alpha_r_dtau_ddelta = 6.125*delta**6*n5*tau**(-0.125) - 9.5*delta**5*n11*tau**3.75*e(-delta**2) -\
                    2.125*delta**5*n8*tau**1.125*e(-delta) - 37.5*delta**4*n12*tau**11.5*e(-delta**3) +\
                    10.625*delta**4*n8*tau**1.125*e(-delta) + 19.0*delta**3*n11*tau**3.75*e(-delta**2) -\
                    13.0*delta**2*n10*tau**5.5*e(-delta**2) + 0.75*delta**2*n4*tau**(-0.75) -\
                    2*delta**2*n7*tau*e(-delta) - 7.0*delta**2*n9*tau**2.5*e(-delta**2) +\
                    25.0*delta*n12*tau**11.5*e(-delta**3) - 2.375*delta*n6*tau**1.375*e(-delta) +\
                    4*delta*n7*tau*e(-delta) + 0.25*n1*tau**(-0.75) + 6.5*n10*tau**5.5*e(-delta**2) +\
                    1.25*n2*tau**0.25 + 1.5*n3*tau**0.5 + 2.375*n6*tau**1.375*e(-delta) + 3.5*n9*tau**2.5*e(-delta**2)


        alpha_0 = a1 + a2 * tau + log(delta) + (c0 - 1) * log(tau) + \
                v1 * log(1 - e((-u1 * tau)/T_c)) + \
                v2 * log(1 - e((-u2 * tau)/T_c)) + \
                v3 * log(1 - e((-u3 * tau)/T_c))

        # Partial alpha(0) with respect to partial delta 
        ##dalpha_0_ddelta = 1/delta 

        # 2nd Partial alpha(0) with respect to partial delta 
        ##d2alpha_0_ddelta2 = -1/delta**2

        # Partial alpha(0) with respect to partial tau
        dalpha_0_dtau = a2 + (c0 - 1)/tau + u1*v1*e(-tau*u1/T_c)/(T_c*(1 - e(-tau*u1/T_c))) +\
                    u2*v2*e(-tau*u2/T_c)/(T_c*(1 - e(-tau*u2/T_c))) + u3*v3*e(-tau*u3/T_c)/(T_c*(1 - e(-tau*u3/T_c))) 

        # 2nd Partial alpha(0) with respect to partial tau
        d2alpha_0_dtau2 = -((c0 - 1)/tau**2 + u1**2*v1*e(-tau*u1/T_c)/(T_c**2*(1 - e(-tau*u1/T_c))) +\
            u1**2*v1*e(-2*tau*u1/T_c)/(T_c**2*(1 - e(-tau*u1/T_c))**2) + u2**2*v2*e(-tau*u2/T_c)/(T_c**2*(1 - e(-tau*u2/T_c))) +\
            u2**2*v2*e(-2*tau*u2/T_c)/(T_c**2*(1 - e(-tau*u2/T_c))**2) + u3**2*v3*e(-tau*u3/T_c)/(T_c**2*(1 - e(-tau*u3/T_c))) +\
            u3**2*v3*e(-2*tau*u3/T_c)/(T_c**2*(1 - e(-tau*u3/T_c))**2))

        dalpha_drho = R*T*\
            (1/rho \
                + n1*tau**0.25/rho_c \
                + n2*tau**1.25/rho_c \
                + n3*tau**1.5/rho_c \
                + 3*n4*rho**2*tau**0.25/rho_c**3 \
                + 7*n5*rho**6*tau**0.875/rho_c**7
                + n6*tau**2.375*e(-delta)/rho_c \
                - n6*rho*tau**2.375*e(-delta)/rho_c**2 \
                + 2*n7*rho*tau**2*e(-delta)/rho_c**2 \
                - n7*rho**2*tau**2*e(-delta)/rho_c**3 \
                + 5*n8*rho**4*tau**2.125*e(-delta)/rho_c**5 \
                - n8*rho**5*tau**2.125*e(-delta)/rho_c**6 \
                + n9*tau**3.5*e(-delta**2)/rho_c \
                - 2*n9*rho**2*tau**3.5*e(-delta**2)/rho_c**3 \
                + n10*tau**6.5*e(-delta**2)/rho_c \
                - 2*n10*rho**2*tau**6.5*e(-delta**2)/rho_c**3 \
                + 4*n11*rho**3*tau**4.75*e(-delta**2)/rho_c**4 \
                - 2*n11*rho**5*tau**4.75*e(-delta**2)/rho_c**6 \
                + 2*n12*rho*tau**12.5*e(-delta**3)/rho_c**2 \
                - 3*n12*rho**4*tau**12.5*e(-delta**3)/rho_c**5 )

        #Compressibility Factor 
        Z = 1 + delta * dalpha_r_ddelta

        #internal energy 
        u = tau * (dalpha_0_dtau + dalpha_r_dtau) * R * T

        #enthalpy
        h = R* T * (tau * (dalpha_0_dtau + dalpha_r_dtau) + delta * dalpha_r_ddelta + 1)

        #entropy
        s = R * (tau * (dalpha_0_dtau + dalpha_r_dtau)-alpha_0-alpha_r)

        #constant volume heat capacity (isochoric)
        C_v = R* (-tau**2 * (d2alpha_0_dtau2 + d2alpha_r_dtau2))

        #constant pressure heat capacity (isobaric)
        C_p = R * (C_v / R + (1 + delta * dalpha_r_ddelta - delta * tau * d2alpha_r_dtau_ddelta)**2/\
            (1 + 2 * delta * dalpha_r_ddelta + (delta**2) * d2alpha_r_ddelta2))

        #pressure
        P = rho**2 * (dalpha_drho)

        return Z, u, h, s, C_v, C_p, P

    def UnSat_Thermal_Cond(self, T, P):
        
        """
        Calculates the termal conductivity of unsaturated N2O between 277K and 444K
        From "Thermal Conductivity of Fluids. Nitrous Oxide"
        By: G. NEAL RICHTER and B. H. SAGE
        """

        K_Table = [\
            [277.594444, [
                [101352.932,	0.0155766120],
                [1378951.458,	0.0168054336],
                [2757902.916,	0.0189169299],
                [4136854.374,	0.1019402721],
                [5515805.832,	0.1036363920],
                [6894757.290,	0.1051940532],
                [10342135.935,	0.1086901373],
                [13789514.580,	0.1118919964],
                [17236893.225,	0.1148688600],
                [20684271.870,	0.1177245723],
                [24131650.515,	0.1205456698],
                [27579029.160,	0.1233667673],
                [31026407.805,	0.1261532501],
                [34473786.450,	0.1290781917] ]\
            ],
            [310.927778, [
                [101352.932,	0.0178957965],
                [1378951.458,	0.0186573197],
                [2757902.916,	0.0200592148],
                [4136854.374,	0.0218245642],
                [5515805.832,	0.0248533499],
                [6894757.290,	0.0334897159],
                [10342135.935,	0.0735562235],
                [13789514.580,	0.0829368054],
                [17236893.225,	0.0884924637],
                [20684271.870,	0.0922654653],
                [24131650.515,	0.0950865628],
                [27579029.160,	0.0975788207],
                [31026407.805,	0.0995691656],
                [34473786.450,	0.1017845059] ]\
            ],
            [344.261111, [
                [101352.932,	0.0206130499],
                [1378951.458,	0.0211841924],
                [2757902.916,	0.0219803303],
                [4136854.374,	0.0236245282],
                [5515805.832,	0.0258225613],
                [6894757.290,	0.0283148192],
                [10342135.935,	0.0356704415],
                [13789514.580,	0.0502086128],
                [17236893.225,	0.0625660583],
                [20684271.870,	0.0696274558],
                [24131650.515,	0.0744042835],
                [27579029.160,	0.0777965234],
                [31026407.805,	0.0802887813],
                [34473786.450,	0.0823829703] ]\
            ],
            [377.594444, [
                [101352.932,	0.0233303034],
                [1378951.458,	0.0237110650],
                [2757902.916,	0.0243687442],
                [4136854.374,	0.0254244923],
                [5515805.832,	0.0266533139],
                [6894757.290,	0.0277782915],
                [10342135.935,	0.0325378118],
                [13789514.580,	0.0389242227],
                [17236893.225,	0.0458817761],
                [20684271.870,	0.0527874074],
                [24131650.515,	0.0592949698],
                [27579029.160,	0.0644698665],
                [31026407.805,	0.0679313358],
                [34473786.450,	0.0712716537] ]\
            ],
            [410.927778,[
                [101352.932,	0.0260129421],
                [1378951.458,	0.0262898596],
                [2757902.916,	0.0269648462],
                [4136854.374,	0.0277090621],
                [5515805.832,	0.0285225074],
                [6894757.290,	0.0293705673],
                [10342135.935,	0.0325205044],
                [13789514.580,	0.0370204146],
                [17236893.225,	0.0418837790],
                [20684271.870,	0.0471798271],
                [24131650.515,	0.0520604989],
                [27579029.160,	0.0561104180],
                [31026407.805,	0.0596930388],
                [34473786.450,	0.0624968289] ]\
            ],
            [444.261111, [
                [101352.932,	0.0286955808],
                [1378951.458,	0.0290071131],
                [2757902.916,	0.0295609482],
                [4136854.374,	0.0301840126],
                [5515805.832,	0.0308589992],
                [6894757.290,	0.0315686004],
                [10342135.935,	0.0337147114],
                [13789514.580,	0.0363281207],
                [17236893.225,	0.0396857460],
                [20684271.870,	0.0433895182],
                [24131650.515,	0.0468509875],
                [27579029.160,	0.0501220760],
                [31026407.805,	0.0532720131],
                [34473786.450,	0.0564219502] ]\
            ]\
        ]

        def K_from_P(temp_idx, P):
            """
            Returns the pressure dependent thermal conductivity of unsaturated N2O at one or two points from data tables. 
            """
            for pressure in range(0, 14):
                if pressure == 13: #   executes if when loop is on the last data set
                    if K_Table[temp_idx][1][pressure][0] == P:
                        P_data = 0
                        K_data = [K_Table[temp_idx][1][pressure][1]]
                        return P_data, K_data
                elif K_Table[temp_idx][1][pressure][0] <= P <= K_Table[temp_idx][1][pressure + 1][0]:
                    P_data = [K_Table[temp_idx][1][pressure][0], K_Table[temp_idx][1][pressure + 1][0]]
                    K_data = [K_Table[temp_idx][1][pressure][1], K_Table[temp_idx][1][pressure + 1][1]]
                    return P_data, K_data

        def K_from_T(T, P):
            """
            Finds the thermal conductivity at current T & P through single or double interpolation. 
            """

            for temp_idx in range(0, 6):
                if temp_idx == 5:    
                    P_data, K_data = K_from_P(temp_idx, P)
                    if P_data == 0: # checks for single K value returned 
                        k = K_data[0]
                        return k
                    else:
                        k = np.interp(P, P_data, K_data)
                        return k
                else:
                    # finds thermal cond. at current T & P for lower temp array 
                    P_data_low, K_data_low = K_from_P(temp_idx, P)    
                    K_low = np.interp(P, P_data_low, K_data_low)
                    # finds thermal cond. at current T & P for higher temp array 
                    P_data_high, K_data_high = K_from_P(temp_idx + 1, P)    
                    K_high = np.interp(P, P_data_high, K_data_high)
                    
                    # creates T & k array for interpolation
                    T_data = [K_Table[temp_idx][0], K_Table[temp_idx + 1][0]]
                    K_data = [K_low, K_high]

                    k = np.interp(T, T_data, K_data)
                    return k        

        if not 277.594444 <= T <= 444.261111:
            raise RuntimeWarning("Input temperature exceedes input limit:\n 277.594444<= T <= 444.261111")
        elif not 101352.932 <= P <= 34473786.450:
            raise RuntimeWarning("Input pressure exceeds input limit:\n 101352.932<= P <= 34473786.450")  
        else:
            print("Finding Sup Crit K")
            k = K_from_T(T, P)

        return k

    def Unsat_Visc(self, T, P):
        """
        Calculates the viscosity of unsaturated N2O between 298K and 398K
        From "Viscosity of Gaseous Nitrous Oxide from 298.15 K to 398.15 K 
                    at Pressures up to 25 MPa"
        By: Takahashi et. al. 
        """

        V_Table = [\
                    [298.15,[ 
                        [102500,	1.490042E-05],
                        [105000,	1.490128E-05],
                        [115000,	1.490474E-05],
                        [130000,	1.490993E-05],
                        [145000,	1.491512E-05],
                        [160000,	1.491964E-05],
                        [190000,	1.490766E-05],
                        [210000,	1.490006E-05],
                        [250000,	1.490291E-05],
                        [300000,	1.490647E-05],
                        [350000,	1.491002E-05],
                        [400000,	1.491301E-05],
                        [500000,	1.491900E-05],
                        [600000,	1.492972E-05],
                        [700000,	1.494182E-05],
                        [800000,	1.495700E-05],
                        [900000,	1.497292E-05],
                        [1000000,	1.499318E-05],
                        [1250000,	1.503012E-05],
                        [1500000,	1.506667E-05],
                        [1750000,	1.513004E-05],
                        [2000000,	1.518741E-05],
                        [2250000,	1.524793E-05],
                        [2500000,	1.531221E-05],
                        [2750000,	1.538563E-05],
                        [3000000,	1.548272E-05],
                        [3250000,	1.560132E-05],
                        [3500000,	1.573206E-05],
                        [3750000,	1.588373E-05],
                        [4000000,	1.602321E-05],
                        [4250000,	1.623704E-05],
                        [4500000,	1.645714E-05],
                        [4750000,	1.668757E-05],
                        [5000000,	1.699236E-05],
                        [5250000,	1.735528E-05],
                        [5500000,	1.785124E-05] ] ],
                    
                    
                    [323.15, [
                        [102500,	1.601000E-05],
                        [105000,	1.601000E-05],
                        [115000,	1.601000E-05],
                        [130000,	1.601000E-05],
                        [145000,	1.601000E-05],
                        [160000,	1.601000E-05],
                        [190000,	1.601000E-05],
                        [210000,	1.601070E-05],
                        [250000,	1.601596E-05],
                        [300000,	1.602254E-05],
                        [350000,	1.602912E-05],
                        [400000,	1.603570E-05],
                        [500000,	1.604886E-05],
                        [600000,	1.606122E-05],
                        [700000,	1.606920E-05],
                        [800000,	1.607718E-05],
                        [900000,	1.609541E-05],
                        [1000000,	1.611927E-05],
                        [1250000,	1.614395E-05],
                        [1500000,	1.619490E-05],
                        [1750000,	1.625650E-05],
                        [2000000,	1.631419E-05],
                        [2250000,	1.635030E-05],
                        [2500000,	1.642392E-05],
                        [2750000,	1.650486E-05],
                        [3000000,	1.658965E-05],
                        [3250000,	1.666991E-05],
                        [3500000,	1.677364E-05],
                        [3750000,	1.688620E-05],
                        [4000000,	1.699756E-05],
                        [4250000,	1.710846E-05],
                        [4500000,	1.722789E-05],
                        [4750000,	1.738240E-05],
                        [5000000,	1.754097E-05],
                        [5250000,	1.771731E-05],
                        [5500000,	1.791126E-05] ] ],
                    
                    [348.15, [
                        [102500,	1.715047E-05],
                        [105000,	1.715164E-05],
                        [115000,	1.715633E-05],
                        [130000,	1.716336E-05],
                        [145000,	1.717038E-05],
                        [160000,	1.717741E-05],
                        [190000,	1.719147E-05],
                        [210000,	1.720000E-05],
                        [250000,	1.720000E-05],
                        [300000,	1.720000E-05],
                        [350000,	1.719983E-05],
                        [400000,	1.719643E-05],
                        [500000,	1.719134E-05],
                        [600000,	1.721578E-05],
                        [700000,	1.724013E-05],
                        [800000,	1.725468E-05],
                        [900000,	1.726923E-05],
                        [1000000,	1.728374E-05],
                        [1250000,	1.732908E-05],
                        [1500000,	1.738561E-05],
                        [1750000,	1.742317E-05],
                        [2000000,	1.745061E-05],
                        [2250000,	1.752038E-05],
                        [2500000,	1.759420E-05],
                        [2750000,	1.760045E-05],
                        [3000000,	1.771306E-05],
                        [3250000,	1.780704E-05],
                        [3500000,	1.789357E-05],
                        [3750000,	1.798295E-05],
                        [4000000,	1.807526E-05],
                        [4250000,	1.817593E-05],
                        [4500000,	1.829000E-05],
                        [4750000,	1.840070E-05],
                        [5000000,	1.849901E-05],
                        [5250000,	1.863807E-05],
                        [5500000,	1.876520E-05] ] ],
                    
                    [373.15, [
                        [102500,	1.835005E-05],
                        [105000,	1.835021E-05],
                        [115000,	1.835086E-05],
                        [130000,	1.835183E-05],
                        [145000,	1.835280E-05],
                        [160000,	1.835377E-05],
                        [190000,	1.835571E-05],
                        [210000,	1.835700E-05],
                        [250000,	1.835959E-05],
                        [300000,	1.836531E-05],
                        [350000,	1.837141E-05],
                        [400000,	1.837750E-05],
                        [500000,	1.838968E-05],
                        [600000,	1.840126E-05],
                        [700000,	1.841283E-05],
                        [800000,	1.842580E-05],
                        [900000,	1.844107E-05],
                        [1000000,	1.845634E-05],
                        [1250000,	1.850291E-05],
                        [1500000,	1.855516E-05],
                        [1750000,	1.860037E-05],
                        [2000000,	1.864941E-05],
                        [2250000,	1.870914E-05],
                        [2500000,	1.876390E-05],
                        [2750000,	1.884613E-05],
                        [3000000,	1.889556E-05],
                        [3250000,	1.896649E-05],
                        [3500000,	1.904194E-05],
                        [3750000,	1.910063E-05],
                        [4000000,	1.917950E-05],
                        [4250000,	1.927255E-05],
                        [4500000,	1.936622E-05],
                        [4750000,	1.945400E-05],
                        [5000000,	1.954040E-05],
                        [5250000,	1.965858E-05],
                        [5500000,	1.975561E-05] ] ],
                    
                    [398.15, [
                        [102500,	1.833997E-05],
                        [105000,	1.833955E-05],
                        [115000,	1.833789E-05],
                        [130000,	1.833539E-05],
                        [145000,	1.833289E-05],
                        [160000,	1.833039E-05],
                        [190000,	1.832540E-05],
                        [210000,	1.832206E-05],
                        [250000,	1.831540E-05],
                        [300000,	1.830708E-05],
                        [350000,	1.830148E-05],
                        [400000,	1.831133E-05],
                        [500000,	1.833103E-05],
                        [600000,	1.835011E-05],
                        [700000,	1.835318E-05],
                        [800000,	1.835625E-05],
                        [900000,	1.835932E-05],
                        [1000000,	1.836594E-05],
                        [1250000,	1.838497E-05],
                        [1500000,	1.840387E-05],
                        [1750000,	1.842296E-05],
                        [2000000,	1.844354E-05],
                        [2250000,	1.846661E-05],
                        [2500000,	1.849965E-05],
                        [2750000,	1.853012E-05],
                        [3000000,	1.855649E-05],
                        [3250000,	1.858220E-05],
                        [3500000,	1.860257E-05],
                        [3750000,	1.862576E-05],
                        [4000000,	1.866576E-05],
                        [4250000,	1.870421E-05],
                        [4500000,	1.873341E-05],
                        [4750000,	1.876972E-05],
                        [5000000,	1.881535E-05],
                        [5250000,	1.885563E-05],
                        [5500000,	1.888961E-05] ] ] ]

        def V_from_P(temp_idx, P):
            """
            Returns the pressure dependent viscosity of unsaturated N2O at one or two points from data tables. 
            """
            print("V from P attempt")
            for pressure in range(0, 36):
                if pressure == 13: #   executes if when loop is on the last data set
                    if V_Table[temp_idx][1][pressure][0] == P:
                        print("V from P if attempt")
                        P_data = 0
                        V_data = [V_Table[temp_idx][1][pressure][1]]
                        print("V from P Found")  
                        return P_data, V_data
                elif V_Table[temp_idx][1][pressure][0] <= P <= V_Table[temp_idx][1][pressure + 1][0]:
                    print("V from P attempt @ else")
                    P_data = [V_Table[temp_idx][1][pressure][0], V_Table[temp_idx][1][pressure + 1][0]]
                    V_data = [V_Table[temp_idx][1][pressure][1], V_Table[temp_idx][1][pressure + 1][1]]
                    print("V from P Found")  
                    return P_data, V_data
  
        def V_from_T(T, P):
            """
            Finds the Viscosity at current T & P through single or double interpolation. 
            """
            for temp_idx in range(0, 5):
                if temp_idx == 4:    
                    P_data, V_data = V_from_P(temp_idx, P)
                    if P_data == 0: # checks for single K value returned 
                        V = V_data[0]
                        print("V from T Found")  
                        return V
                    else:
                        V = np.interp(P, P_data, V_data)
                        print("V from T Found")  
                        return V
                else:
                    # finds thermal cond. at current T & P for lower temp array 
                    print("V from T attempt @ else")
                    P_data_low, V_data_low = V_from_P(temp_idx, P)    
                    print("L1012")
                    V_low = np.interp(P, P_data_low, V_data_low)
                    print("L1014")
                    # finds thermal cond. at current T & P for higher temp array 
                    P_data_high, V_data_high = V_from_P(temp_idx + 1, P)    
                    print("L1017")
                    V_high = np.interp(P, P_data_high, V_data_high)
                    print("L1019")

                    # creates T & k array for interpolation
                    T_data = [V_Table[temp_idx][0], V_Table[temp_idx + 1][0]]
                    V_data = [V_low, V_high]
                    print("L1024")
                    V = np.interp(T, T_data, V_data)
                    print("L1026")
                    print("V from T Found")  
                    return V          

        if not 298.15 <= T <= 398.15:
            raise RuntimeWarning("Input temperature exceedes input limit:\n 298.15 <= T <= 398.15")
        elif not 102500 <= P <= 5500000:
            raise RuntimeWarning("Input pressure exceeds input limit:\n 102500 <= P <= 5500000")  
        else:
            print("V from T Enter")  
            V = V_from_T(T, P)

        return V  

    def EQM_Tank_Solve(self, IC):
        """
        Method where the N2O state equations are solved.
        Assumes heat transfer on the tank wall is instant.
        """
        T_N2O_1 = IC[0]
        U_N2O_1 = IC[1]
        T_wall = IC[2]
        N2O_mass_1 = IC[3]
        P_N2O_1 = IC[4]
        X = IC[5]
        P_out = IC[6]

        # Exception handling for determining N2O properties, skips to next solving iteration if error is encountered
        try:
            # Modified by JS
            if 0 <= X <= 1:
                h1 = PropsSI('H', 'P', P_N2O_1, 'Q', X, 'N2O')
                rho1 = PropsSI('D', 'P', P_N2O_1, 'Q', X, 'N2O')
            else:
                h1 =  float(PropsSI('H', 'T', T_N2O_1, 'P', P_N2O_1, 'N2O'))
                rho1 = float(PropsSI('D', 'T', T_N2O_1, 'P', P_N2O_1, 'N2O'))
            # State 2 is defined at the P of the exit/entry, not the T...JS
            h2 = float(PropsSI('H', 'P', P_out, 'Q', 0, 'N2O'))
            # rho2 = float(PropsSI('D', 'P', P_out, 'Q', 0, 'N2O'))    # JS 3/11/2011
            T2 = float(PropsSI('T', 'P', P_out, 'Q', 0, 'N2O'))    # JS 3/11/2011
            print("N2O Prop Find Attempt")
            _, rho2, mu_n2o_2, _, _ = self.N2O_prop(T2, 0)    # JS 3/11/2011
            print("h2, T2 found | L1044")

        except Exception as excpt:
            print(excpt)  # Added by JS
            raise
        # Calls Mass Flux method which returns mass flow rate.
        # Negative mass flow indicates flow out of the tank.
        # m_dot = self.Mass_flux(A_out, P_N2O, h1, rho1, P_out, h2, rho2)
        # JS method finds mdot after viscosity, below

        T_film_wall_outside = 0.5 * (T_wall + self.T_atm) # K   | Film temperature between tank wall and atmosphere
        T_film_wall_inside = 0.5 * (T_wall + T_N2O_1)  # K   | Film temperature between tank wall and saturated N2O

        # Calculates and unpacks. Simplified by JS
        Cp_air, rho_air, mu_air, k_air, beta_air = self.Air_prop(T_film_wall_outside)
        
        if 182.33 < T_film_wall_inside < 309.52:
            if 0 <= X <= 1:
                Cp_n2o, N2O_density, mu_n2o, k_n2o, beta_n2o = self.N2O_prop(T_film_wall_inside, X)
            else:
                Cp_n2o, N2O_density, mu_n2o, k_n2o, beta_n2o = self.N2O_Unsat(T_film_wall_inside, P_N2O_1, N2O_mass_1)
        else:
                Cp_n2o, N2O_density, mu_n2o, k_n2o, beta_n2o = self.N2O_Unsat(T_film_wall_inside, P_N2O_1, N2O_mass_1)

        # JS modification - use new method for filling, old method for draining
        if self.Fill:
            m_dot = self.Mass_flow_JS([P_N2O_1, P_out], rho2, mu_n2o_2)    # JS 3/11/2011
            print("M_dot calculated | L1073")
        else:
            m_dot = self.Mass_flux(self.A_tank_exit, P_N2O_1, h1, rho1, P_out, h2, rho2)

        ######
        ######

        print("rho = {}".format(N2O_density))
        print("Cp = {}".format(Cp_n2o))
        print("beta = {}".format(beta_n2o))

        dQ_Ra_nominator = Cp_n2o * N2O_density**2 * g * beta_n2o * abs(T_wall - T_N2O_1) * self.L_tank**3
        dQ_Ra_denominaor = mu_n2o * k_n2o
        dQ_Ra = dQ_Ra_nominator / dQ_Ra_denominaor
        if dQ_Ra < 0:
            dQ_Ra = -1 * dQ_Ra
        
        dQ_Nu = 0.21 * dQ_Ra**(2/5)

        dQ_wall_inside_dt = dQ_Nu * (k_n2o/self.L_tank) * self.A_tank_in * (T_wall - T_N2O_1)

        ######
        ######

        dQ_wall_outside_dt = ((0.59 * (  (Cp_air * (rho_air )**2 * g * beta_air * abs(T_wall - self.T_atm) * self.L_tank**3)\
               / (mu_air * k_air ) )**(1/4) * (k_air/self.L_tank)) * self.A_tank_out * (T_wall - self.T_atm))

        ######
        ######

        dT_wall_dt =  (dQ_wall_outside_dt - dQ_wall_inside_dt) /(self.m_wall * C_wall)
        dT_wall = dT_wall_dt * t_step

        ######
        ######

        dQ_in_dt = dQ_wall_inside_dt #+ dQ_wall_cond_dt
        dQ = dQ_in_dt * t_step

        ######
        ######

        dm_dt = m_dot # This is dm/dt......JS
        dm = dm_dt * t_step # This is dm....JS

        ###############

        try:
            dU_tot_dt = m_dot * h2 + dQ_in_dt  # This should use the enthalpy at the inlet/exit, h2.........JS
            dU_tot = dU_tot_dt * t_step
            print("DU TOT dt = {}".format(dU_tot_dt))
            print("DU TOT = {}".format(dU_tot))
        except:
            return 'next'

        ###############
        #change of density with respect to time
        ###############

        U_N2O_2 = U_N2O_1 + dU_tot     # U at state 2
        N2O_mass_2 = N2O_mass_1 + dm  # m at state 2
        T_wall_2 = T_wall + dT_wall

        ###############
        ###############


        print("T_N2O_1 = {}".format(T_N2O_1))

        if N2O_mass_1 == 0:
            T_N2O_2 = T_N2O_1
            X_2 = 2
            P_N2O_2 = float(PropsSI('P', 'T', T_N2O_2, 'D', dm/self.V_tank_inside, 'N2O')) + self.P_atm

        elif T_N2O_1 < 309.52:
            temp_array = self.Quality(U_N2O_2, N2O_mass_2)

            if temp_array == 'next':
                print('Exit @ EQM_Tank_Solve "Try" | L1162')
                return temp_array

            elif temp_array == True:
                print('Exit @ EQM_Tank_Solve "self.end_loop" | L1166')
                self.end_loop = True #exits loop
                return self.end_loop

            T_N2O_2, X_2, P_N2O_2 = temp_array

        elif T_N2O_1 >= 309.52:

            print("N2O Phase")
            print(CP.PhaseSI('T',T_N2O_1,'D',N2O_mass_2/self.V_tank_inside,'N2O'))

            T_N2O_2 = self.Find_T_from_u(U_N2O_2, N2O_mass_2/self.V_tank_inside)

            _, _, _, _, _, _, P_N2O_2 = self.Helmholtz_alpha_eqs(T_N2O_2, (N2O_mass_2/self.V_tank_inside))
                        
            #sets the quality outside the range of (0,1) to indicate a non-saturated state
            if T_N2O_2 > 309.52:
                X_2 = 2
            elif T_N2O_2 < 182.33:
                X_2 = -1

        V_out = abs(m_dot)/(rho2 * self.A_tank_exit)  # JS version

        print("\nNew Pressure = {:.2f} | L1198\n".format(P_N2O_2))
        return T_N2O_2, U_N2O_2, T_wall_2, N2O_mass_2, P_N2O_2, X_2

    def Quality(self, U, mass):
        """
        Determins the Quality, Temperature, & Pressure of the N2O tank.
        Attempts to use ideal gas law if the density << theoretical density.
        Calculates compressability factor for find Pressure in the tank if this
        is the case.
        Otherwise, uses internal energy find quality and volume constrain to find
        density and pressure.
        """
        
        Volume = self.V_tank_inside

        def Quality_JS(T, V, U, m):    # JS 3/11/2021
            """
            Uses rootfinding to find a temperature such that it satisfies the following:
                1) V/m = (1-x)/rho_l + x/rho_v
                and
                2) X = (U/m- u_l)/(u_v - u_l)
            Created by JS on 3/11/2021
            Enthalpy is typically calculated relative to some reference point.
            I'm guessing you were getting problems because you were mixing
            CoolProps's enthalpy with some other polynomial fit.
            """
            print("Finding Quality_JS | L1259")

            if T <= 182.33:
                print ("Quality Temp = {}".format(T))
                print("T < 182.33: Returned X = 1000 | L1261")
                return 1000
            elif T >= 309.52:
                print ("Quality Temp = {}".format(T))
                print("T > 309.52: Returned X = 0 | L1264")
                return 1
            else:
                print ("Quality Temp = {}".format(T))
                u_l = PropsSI('U', 'T', T, 'Q', 0, 'N2O')
                u_v = PropsSI('U', 'T', T, 'Q', 1, 'N2O')

                x = (U/m - u_l)/(u_v - u_l)
                x2 = quality_from_T_rho(T[0], m/V)

                print("Quality_JS Found | L1278")
                return x - x2


        def quality_from_T_rho(T, rho): # JS 3/11/2021.
            print('Finding quality_from_T_rho | L1283')
            rho_l = PropsSI('D', 'T', T, 'Q', 0, 'N2O')
            rho_v = PropsSI('D', 'T', T, 'Q', 1, 'N2O')
            x = (1/rho - 1/rho_l)/(1/rho_v - 1/rho_l)
            print("Calculated Quality = {} | L1289".format(x))
            return x

        args = (Volume, U, mass)
        eq_sol = root(Quality_JS, 260 , args=args)  # JS Note: This assumes the tank is a saturated mixture.  # JS 3/11/2021
        print(eq_sol.message)
        print(eq_sol.x)
        Temp_sol = eq_sol.x[0]

        if 182.33 < Temp_sol < 309.52:
            # X_sol = PropsSI("Q", "T", Temp_sol, "D", mass/Volume, "N2O")  # JS 3/11/2021
            X_sol = quality_from_T_rho(Temp_sol, mass/Volume)  # JS 3/11/2021
            print(X_sol)
            if 0 <= X_sol <= 1:
                P_N2O = PropsSI('P', 'T', Temp_sol, 'Q', X_sol, 'N2O')
            else:
                P_N2O = PropsSI('P', 'T', Temp_sol, 'D', mass/Volume, 'N2O')
        elif Temp_sol >= 309.52:
            X_sol = 2
            _, _, _, _, _, _, P_N2O = self.Helmholtz_alpha_eqs(Temp_sol, (mass/Volume))

        ###### Block below # JS 3/11/2021

        print("Quality = {:.4f} | L1321".format(X_sol))
        print("Calculated P = {:.2f} | L1322".format(P_N2O))

        return Temp_sol, X_sol, P_N2O

    def Pipe_Solve(self, mdot, T, X, P_in, Vel_):
        """
        Calculates the pressure at the end of the pipe, assumes constant velocity along the stream line.
        Minor loss coefficient found in Pipe_Info() method.
        """
        # Note::Assumes constant velocity throughout pipe
        # Accuracy can be improved
        # Assumes Saturated liquid flow inside pipe
        # Assumes constant height, will need to be changed for accuracy
        mdot = abs(mdot)
        N2O_Coeff = self.N2O_prop(T, 0) # I changed this to use sat. liquid......JS
        mu = N2O_Coeff[2]
        rho = PropsSI('D', 'T', T, 'Q', 0, 'N2O') #kg/m^3
        Vel_in = mdot/(rho * self.A_pipe) #m/s

        Re = (rho * Vel_in * self.D_pipe )/ mu

        # The last argument, Vel_, isn't used and can be removed.......JS
        print("Re = {}".format(Re))
        if Re <= 2000:
            f = 64/Re
        else:
            # Haaland Equation
            f = ( 1 / \
                    (-1.8*log( ( (self.rough/self.D_pipe)/(3.7) )**1.11 + \
                        6.9/Re)))**2
        # Major Loss in Pipe
        H_major = f * self.L_pipe/self.D_pipe * Vel_in**2/(2*g)
        print("H_major = {}".format(H_major))

        # Minor Loss in Pipe
        H_minor = self.K_pipe * Vel_in**2/(2*g)
        # Total Loss
        print("H_minor = {}".format(H_minor))
        H_loss = H_major + H_minor
        print("H_loss = {}".format(H_loss))
        P_out = P_in - H_loss * rho * g

        return P_out

    def Visualize(self):
        """
        Info:
        Method plots the N2O state data.
        """
        print("Program terminated")

        fig, ax = plt.subplots(figsize=(8,8))
        plt.ylabel('Temperature (K)')
        plt.xlabel('Time (s)')
        line1, = ax.plot(self.t_plt, self.plot_0, label='Saturated N2O Temperature')
        #line2, = ax.plot(self.t_plt, self.plot_2, label='Tank Wall Temperature')
        plt.grid(True)
        plt.title('Temperature Data')
        plt.legend()
        plt.show()

        fig, ax = plt.subplots(figsize=(8,8))
        plt.ylabel('Internal Energy (J)')
        plt.xlabel('Time (s)')
        line1, = ax.plot(self.t_plt, self.plot_1, label='N2O Internal Energy')
        plt.grid(True)
        plt.title('Saturated N2O Internal Energy')
        plt.legend()
        plt.show()

        fig, ax = plt.subplots(figsize=(8,8))
        plt.ylabel('N2O Mass (kg)')
        plt.xlabel('Time (s)')
        line1, = ax.plot(self.t_plt, self.plot_3, label='N2O Mass')
        plt.grid(True)
        plt.title('N2O Tank Mass')
        plt.legend()
        plt.show()

        fig, ax = plt.subplots(figsize=(8,8))
        plt.ylabel('Pressure (kPa)')  # JS 3/11/2021
        plt.xlabel('Time (s)')
        line1, = ax.plot(self.t_plt, [x/1000 for x in self.plot_4], label='Pressure') # JS 3/11/2021
        plt.grid(True)
        plt.title('N2O Tank Pressure')
        plt.legend()
        plt.show()

        fig, ax = plt.subplots(figsize=(8,8))
        plt.ylabel('Quality')
        plt.xlabel('Time (s)')
        line1, = ax.plot(self.t_plt, self.plot_5, label='Quality')
        plt.ylim(0,1)
        plt.grid(True)
        plt.title('N2O Quality')
        plt.legend()
        plt.show()

        # fig, ax = plt.subplots(figsize=(8,8))
        # plt.ylabel('Mass Flow Rate (kg/s)')
        # plt.xlabel('Time (s)')
        # line1, = ax.plot(self.t_plt, self.plot_7, label='Mass Flow Rate')
        # plt.grid(True)
        # plt.title('Tank Exit N2O Mass Flow Rate')
        # plt.legend()
        # plt.show()

        # fig, ax = plt.subplots(figsize=(8,8))
        # plt.ylabel('Pressure (Pa)')
        # plt.xlabel('Time (s)')
        # line1, = ax.plot(self.t_plt, self.plot_6, label='Exit Pressure')
        # plt.grid(True)
        # plt.title('N2O Exit Pressure')
        # plt.legend()
        # plt.show()

    def Find_T_from_u(self, U, rho):
        """
        Uses the finding with the Helmholtz functional equations to find the 
        Temperature at a given internal energy and density
        From: Short Fundamental Equations of State for 20 Industrial Fluids
        By: E.W. Lemmon & R. Span
        """
        print("Find_T_from_u entered | L1447")
        def T_Solve(T, U, rho):

            T_c = 309.52
            rho_c = 452
            R_m  = 8.31446261815324 # J/(mol K) Molar Gas Constant
            m_w = (44.0128/1000) #kg/mol N2O Molecular Weight 
            R = R_m/m_w #J/(kg K) N2O Specific Gas Constant

            delta = rho/rho_c
            tau = T / T_c

            #Alpha_r constants
            n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12 = \
                0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032, \
                    -0.0036985, -0.23263, -0.00042859, -0.042810, -0.023038

            #Alpha_0 constants
            v1, u1, v2, u2, v3, u3 = 2.1769, 879.0, 1.6145, 2372.0, 0.48393, 5447.0
            c0, a1, a2 = 3.5, -4.4262736272, 4.3120475243

            alpha_r = n1 * delta * (tau**0.25) +\
                        n2 * delta * (tau**1.25) +\
                        n3 * delta * (tau**1.5) + \
                        n4 * (delta**3) * (tau**0.25) + \
                        n5 * (delta**7) * (tau**0.875) + \
                        n6 * delta * (tau**2.375) * (e(-delta)) + \
                        n7 * (delta**2) * (tau**2) * (e(-delta)) + \
                        n8 * (delta**5) * (tau**2.125) * (e(-delta)) + \
                        n9 * (delta) * (tau**3.5) * (e(-delta**2)) + \
                        n10 * delta * (tau**6.5) * (e(-delta**2)) + \
                        n11 * (delta**4) * (tau**4.75) * e(-delta**2) + \
                        n12 * (delta**2) * (tau**12.5) * e(-delta**3)

            # Partial alpha(r) with respect to partial tau
            dalpha_r_dtau = 0.875*delta**7*n5*tau**(-0.125) + 2.125*delta**5*n8*tau**1.125*e(-delta) +\
                        4.75*delta**4*n11*tau**3.75*e(-delta**2) + 0.25*delta**3*n4*tau**(-0.75) +\
                        12.5*delta**2*n12*tau**11.5*e(-delta**3) + 2*delta**2*n7*tau*e(-delta) +\
                        0.25*delta*n1*tau**(-0.75) + 6.5*delta*n10*tau**5.5*e(-delta**2) + 1.25*delta*n2*tau**0.25 +\
                        1.5*delta*n3*tau**0.5 + 2.375*delta*n6*tau**1.375*e(-delta) + 3.5*delta*n9*tau**2.5*e(-delta**2) 

            # Partial alpha(0) with respect to partial tau
            dalpha_0_dtau = a2 + (c0 - 1)/tau + u1*v1*e(-tau*u1/T_c)/(T_c*(1 - e(-tau*u1/T_c))) +\
                        u2*v2*e(-tau*u2/T_c)/(T_c*(1 - e(-tau*u2/T_c))) + u3*v3*e(-tau*u3/T_c)/(T_c*(1 - e(-tau*u3/T_c))) 

            return (tau * (dalpha_0_dtau + dalpha_r_dtau) ) - U/(R * T)

        args = (U, rho)
        Teq_sol = root(T_Solve, 270 , args=args)  # JS Note: This assumes the tank is a saturated mixture.  # JS 3/11/2021
        print(Teq_sol.message)
        print(Teq_sol.x)
        T_sol = Teq_sol.x[0]

        return T_sol


Main()