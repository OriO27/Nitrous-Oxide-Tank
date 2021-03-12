"""
OOrtiz
V: 3.2.2

ADDED:
Comments on the methods to make it easier to follow

TO DO:
Need to fix the pressure gradient at the pipe exit
Need to track total mass transfer through piping
"""

import numpy as np
import matplotlib.pyplot as plt
from math import exp as e
from math import log10 as log
from scipy.optimize import root  # JS added
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
        t_step = 0.1 #s

        k = 10000
        self.N2O_total_mass = 10 #kg  # TEMP CHANGE TO 1 KG to calculate faster..........JS
        self.t = 0.0
        self.P_atm = float(101325) #Pa | atmospheric pressure
        self.Fill = True #True for Fill | False for Empty

        self.end_loop = False
        self.Tank_Geometry()
        self.Pipe_Info()
        self.z = self.Tank_initial(self.Fill)

        self.t_plt = []     # Time
        self.plot_0 = []    # N2O Temp
        self.plot_1 = []    # U_N2O
        self.plot_2 = []    # Wall Temp
        self.plot_3 = []    # N2O Mass
        self.plot_4 = []    # Tank Pressure
        self.plot_5 = []    # Quality
        self.plot_6 = []    # P_out
        #self.plot_7 = []    # P_out

        self.Model(k)
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

            print(solver) #counter check
            EQ_Solver = self.EQM_Tank_Solve(self.z)
            self.t += t_step
            if self.end_loop == True:
                print("Loop break @ solver start")
                break
            elif EQ_Solver == 'next':
                continue
            # PipeSolver = self.Pipe_Solve(EQ_Solver[7], EQ_Solver[0], EQ_Solver[5], EQ_Solver[8] , EQ_Solver[9])  # Commented out by JS, unnecessary with my m_dot method

    # EQ_Solver[8] is P_out, which is the pressure in the tubing into/out of the tank
    # PipeSolver is this pressure minus the pressure drop in the tubing

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
                        
            #EQ_Solver[6]  # Changed by JS
            #self.z[8]  # Changed by JS. The pressure far from the rocket tank is constant - either 1 atm for draining or P of storage tank for filling.
                    
            self.t_plt.append(self.t)
            self.plot_0.append(self.z[0]) # N2O Temp
            self.plot_1.append(self.z[1]) # U_N2O
            self.plot_2.append(self.z[2]) # Wall Temp
            self.plot_3.append(self.z[3]) # N2O Mass
            self.plot_4.append(self.z[4]) # Tank Pressure
            self.plot_5.append(self.z[5]) # Quality
            self.plot_6.append(self.z[6]) # P_out...
            
            #self.plot_8.append(self.z[8]) # P_out
            #print(self.z[3])
            
            # Ends the iteration if the mass falls below 0kg or it goes above 0.5% of the initial tank mass
            if self.z[3] <= 0 or self.z[3] >= 1.005*self.N2O_total_mass:
                print("Loop Break\nMass out of range")
                break

    def Tank_initial(self, Fill):
        """
        Sets the initial temperature conditions for the tank and atmosphere

        Assumes atmospheric pressure is 101kPa

        Solves for the initial internal energy of the tank

        Sets the initial quality and tank pressures -- depends on selected filling or emptying mode
        """

        self.T_N2O = 30 + 273.15 #K | wall temperature at liquid region
        self.T_atm = 279.52 #K
        self.T_wall = 0.5 * (self.T_atm + self.T_N2O)

        if Fill == False: #Sets Initial Mass and Quality if tank is in Emptying Mode
            self.N2O_init_mass = self.N2O_total_mass #kg
            X = 0
            P_Tank = float(PropsSI('P', 'Q', X, 'T', self.T_N2O, 'N2O'))
            P_out = self.P_atm #Pa
        elif Fill == True: #Initial Mass and Quality if tank is in Filling Mode
            self.N2O_init_mass = 0.001 #kg
            X = 1
            P_Tank = self.P_atm #Pa | if Filling Tank Assumed to be at P_atm
            P_out = float(PropsSI('P', 'Q', X, 'T', self.T_N2O, 'N2O'))

        U_init = float(PropsSI('U', 'Q', X, 'T', self.T_N2O, 'N2O') * self.N2O_init_mass )

        return self.T_N2O, U_init, self.T_wall, self.N2O_init_mass, P_Tank, X, P_out

    def Tank_Geometry(self):
        """
        Sets tank geometry, calls Tank_Material() method which enables density, heat capacity, and thermal conductivity of
            selected material. Made a separate method to eventaully allow for the slection of other materials.
        The tank exit diameter is used for the mass flow rate into/outof the tank and outof/into the pipe.
        """

        self.Tank_Material()
        self.L_tank = 1.27 # m (50in)   | Tank Length
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
        self.L_pipe = 3.048 #m | Total length of pipe 10ft
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

        if T >= 309.52:
            print("Input Temperature is ABOVE Critical Temperature!")
            print("{:.2f}K will be set to 309.50K!".format(T))
            T = 309.50
        T_r = T/309.52 #reduced Temp

        if X >= 1:
            X = 1
        elif X <= 0:
            X = 0

        Cp = float(PropsSI('C', 'T', T, 'Q', X, 'N2O'))  # J/kg/K      | Heat Capacity
        rho = float(PropsSI('D', 'T', T, 'Q', X, 'N2O')) # kg/m^3      | Density
        beta = float(PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', T, 'Q', X, 'N2O')) # 1/K   | Expansion Coefficient

        # Dynamic Viscosity
            # EQ 4.9 | kg/(m s)
        theta = (309.52 - 5.24)/(T - 5.24)
        mu_l = (0.0293423 * e(1.6089 * (theta - 1)**(1/3) + 2.0439*(theta - 1)**(4/3)) ) * \
                10**(-6)
            # EQ 4.10 | kg/(m s)
        mu_v = e(3.3281 + -1.18237*(1/T_r - 1)**(1/3) + -0.055155*(1/T_r - 1)**(4/3)) * 10**(-3)

        # Thermal Conductivity
            # EQ 4.11 | W/mK
        k_l = (72.35 * (1 + 1.5*(1-T_r)**(1/3) + -3.5*(1-T_r)**(2/3) + 4.5*(1-T_r))) * 10**(-3)
            # EQ 4.12 | W/mK
        k_v = e(\
                -7.08870 + -0.276962*(1-T_r)**(-2/3) + 2.88672*(1-T_r)**(-1/3) + \
                16.6116*(1 - T_r)**(1/3) + -11.8221*(1 - T_r)**(2/3)) * 10**(-3)

        # JS deleted some unnecessary "if" statements
        mu = float(X * (mu_v - mu_l) + mu_l)
        k = float(X * (k_v - k_l) + k_l)

        return Cp, rho, mu, k, beta

    def Mass_flow_JS(self, Ps, rho, mu):
        """Calculate the mass flow rate

        Governing Equation:
        abs(P1 - P2) = rho * g * head_loss

        Parameters
        ----------
        Ps : tuple of length 2
            Two floats indicating the pressure in the tank and at the far end
            of the tube. Order doesn't matter. Units are Pascals
        rho : float
            Density of the liquid in the tube, kg/m^3
        mu : float
            Viscosity of the liquid in the tube, kg/(m s)

        Returns
        -------
        mdot : float
            Mass flow rate, kg/s
        """
        h_loss = abs(Ps[0] - Ps[1]) / (rho * g)

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
        if self.Fill:
            sign = 1
        else:
            sign = -1

        return velocity * self.A_pipe * rho * sign

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
                print('Exit @ Mass Flow || P2 <= P_atm')
                return self.end_loop

        # Simplified by JS
        G_SPI = Cd * (2 * rho1 * abs(P1-P2))**0.5  # Units are kg/(m^2 s)
        G_HEM = Cd * rho2 * (2 * abs(h1 - h2))**0.5  # Units are kg/(m^2 s)

        P1_sat = P1  # This isn't always true.......JS
        k = ((P1-P2)/(P1_sat - P2))**(.5)  # With the line above, this is always = 1.......JS
        G = float((k*G_SPI + G_HEM)/(1+k))

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

        # A_out = IC[6]
        # #M_out = IC[7]
        # P_out = IC[8]

        # Checks to see if the N2O pressure is +/- e-4 % of saturation range to avoid CoolProp error
        # Reduces N2O pressure to 0.9998% of input value if it's within the range
        # Pressure = float(PropsSI('P', 'T', T_N2O, 'Q', 0 , 'N2O'))
        # if P_N2O >= 0.9999 * Pressure and P_N2O <= 1.0001 * Pressure:
        #     P_N2O = 0.9998 * P_N2O
        # Commented out by JS. This would make the state superheated, never saturated

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
            rho2 = float(PropsSI('D', 'P', P_out, 'Q', 0, 'N2O'))
        except Exception as e:
            print(e)  # Added by JS
            raise


        # Calls Mass Flux method which returns mass flow rate.
        # Negative mass flow indicates flow out of the tank.
        # m_dot = self.Mass_flux(A_out, P_N2O, h1, rho1, P_out, h2, rho2)
        # JS method finds mdot after viscosity, below

        T_film_wall_outside = 0.5 * (T_wall + self.T_atm) # K   | Film temperature between tank wall and atmosphere
        T_film_wall_inside = 0.5 * (T_wall + T_N2O_1)  # K   | Film temperature between tank wall and saturated N2O

        # Calculates and unpacks. Simplified by JS
        Cp_air, rho_air, mu_air, k_air, beta_air = self.Air_prop(T_film_wall_outside)
        Cp_n2o, rho_n2o, mu_n2o, k_n2o, beta_n2o = self.N2O_prop(T_film_wall_inside, X)

        # JS modification - use new method for filling, old method for draining
        if self.Fill:
            m_dot = self.Mass_flow_JS([P_N2O_1, P_out], rho2, mu_n2o)
        else:
            m_dot = self.Mass_flux(self.A_tank_exit, P_N2O_1, h1, rho1, P_out, h2, rho2)


        ######
        ######
        # I didn't look at these dQ terms......JS
        #changed PropsSI to use quality and temp instead of pressure and temp
        dQ_wall_inside_dt = ((0.021 * ((Cp_n2o * PropsSI('D', 'T', T_film_wall_inside, 'P', P_N2O_1, 'N2O')**2 * g * \
            beta_n2o * abs(T_wall - T_N2O_1) * self.L_tank**3  )\
               / (mu_n2o * k_n2o ) )**(2/5) * (k_n2o/self.L_tank)) * self.A_tank_in * (T_wall - T_N2O_1) )

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
        except:
            return 'next'

        ###############

        U_N2O_2 = U_N2O_1 + dU_tot     # U at state 2
        N2O_mass_2 = N2O_mass_1 + dm  # m at state 2
        T_wall_2 = T_wall + dT_wall

        temp_array = self.Quality(U_N2O_2, N2O_mass_1, self.V_tank_inside, P_N2O_1, T_N2O_1)

        if temp_array == 'next':
            print('Exit @ EQM_Tank_Solve "Try"')
            return temp_array
        elif temp_array == True:
            print('Exit @ EQM_Tank_Solve "self.end_loop"')
            self.end_loop = True #exits loop
            return self.end_loop

        T_N2O_2, X_2, P_N2O_2 = temp_array

        # V_out = abs(dm_dot)/(PropsSI('D', 'T', T, 'Q', 0, 'N2O')*self.A_tank_exit)
        
        V_out = abs(m_dot)/(rho2 * self.A_tank_exit)  # JS version
        # return T, dU_tot, dT_wall, dm_dot, P, X, A_out, dm_dot, P_out, V_out
        # The line above returns "dm_dot" in two places. As I mentioned above,
        # This is actually dm. I think the second one should be m_dot instead, as below
        return T_N2O_2, U_N2O_2, T_wall_2, N2O_mass_2, P_N2O_2, X_2

    def Quality(self, U, mass, Volume, P, T1):
        """
        Determins the Quality, Temperature, & Pressure of the N2O tank.
        Attempts to use ideal gas law if the density << theoretical density.
        Calculates compressability factor for find Pressure in the tank if this
        is the case.
        Otherwise, uses internal energy find quality and volume constrain to find
        density and pressure.
        """

        def Quality_V2(T, V, U, m):
            """
            Uses rootfinding to find a temperature such that it satisfies the following:
                1) V/m = (1-x)/rho_l + x/rho_v
                and
                2) X = (U/m- u_l)/(u_v - u_l)
            """
            T_r = T/309.57 #K | Reduced Temperature 
            rho_c = 452 #kg/m^3 | Critical Density
            P_c = 7251000 #Pa | Critical Pressure
            a1, a2, a3, a4 = 1.72328, -0.83950, 0.51060, -0.10412
            b1, b2, b3, b4, b5 = -1.009, -6.28792, 7.50332, -7.90463, 0.639427
            c1, c2, c3, c4, c5 = -200, 116.043, -917.225, 794.779, -589.587
            d1, d2, d3, d4, d5 = -200, 440.055, -459.701, 434.081, 485.338
            f1, f2, f3, f4 = -6.71893, 1.35966, -1.3779, -4.051

            rho_l = rho_c * e(0.386533333333333333333333 * (1-T_r))
                # a1*(1 - T_r)**(1/3) + a2*(1 - T_r)**(2/3) + \
                #                 a3*(1 - T_r) + a4*(1 - T_r)**(4/3) )

            rho_v = rho_c * e( 6.498755* (1/T_r - 1))
                # b1*(1/T_r - 1)**(1/3) + b2*(1/T_r - 1)**(2/3) + \
                #                 b3*(1/T_r - 1) + b4*(1/T_r - 1)**(4/3) + b5 * (1/T_r - 1)**(5/3) )

            h_l = c1 + c2*(1 - T_r)**(1/3) + c3*(1 - T_r)**(2/3) + c4*(1 - T_r) + c5*(1 - T_r)**(4/3)
            h_v = d1 + d2*(1 - T_r)**(1/3) + d3*(1 - T_r)**(2/3) + d4*(1 - T_r) + d5*(1 - T_r)**(4/3)

            P = P_c * e(1/T_r * ( \
                            f1*(1 - T_r) + f2*(1 - T_r)**(3/2) + f3*(1 - T_r)**(5/2) + f4*(1 - T_r)**5 ) )

            u_l = h_l - P*rho_l**(-1)
            u_v = h_v - P*rho_v**(-1)


            return U/m * (rho_l - rho_v) + u_l *(rho_l + rho_v) - (V/m * rho_l * rho_v - rho_v) * (u_v - u_l)

        args = (Volume, U, mass)
        eq_sol = root(Quality_V2, 280 , args=args)
        Temp_sol = eq_sol.x[0]
        
        X_sol = PropsSI("Q", "T", Temp_sol, "D", mass/Volume, "N2O")
        
        if X_sol < 0:
            X_sol = 0

        print(eq_sol.x)

        print(X_sol)
        # if Temp_sol < 182.23 or Temp_sol > 309.52:
        #     #self.Visualize()
        #     self.end_loop = True
        #     print('Exit @ Quality "self.end_loop"')
        #     return self.end_loop
        T_r = Temp_sol/309.57
        f1, f2, f3, f4 = -6.71893, 1.35966, -1.3779, -4.051
        P_N2O = 7251000 * e(1/T_r * ( \
                            f1*(1 - T_r) + f2*(1 - T_r)**(3/2) + f3*(1 - T_r)**(5/2) + f4*(1 - T_r)**5 ) )

        if P_N2O < 101340:
            P_N2O = 101340
        elif P_N2O > 7.245e+06:
            P_N2O = 7.245e+06

        # print("Temp Solution = {:.3f}".format(Temp_sol)) JS
        print("Quality = {:.4f}".format(X_sol))
        print("Calculated P = {}".format(P_N2O))
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
        plt.ylabel('Pressure (Pa)')
        plt.xlabel('Time (s)')
        line1, = ax.plot(self.t_plt, self.plot_4, label='Pressure')
        plt.grid(True)
        plt.title('N2O Tank Pressure')
        plt.legend()
        plt.show()

        fig, ax = plt.subplots(figsize=(8,8))
        plt.ylabel('Quality')
        plt.xlabel('Time (s)')
        line1, = ax.plot(self.t_plt, self.plot_5, label='Quality')
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

        fig, ax = plt.subplots(figsize=(8,8))
        plt.ylabel('Pressure (Pa)')
        plt.xlabel('Time (s)')
        line1, = ax.plot(self.t_plt, self.plot_6, label='Exit Pressure')
        plt.grid(True)
        plt.title('N2O Exit Pressure')
        plt.legend()
        plt.show()

Main()