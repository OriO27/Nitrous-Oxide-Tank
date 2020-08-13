"""
OOrtiz 
N2O Pressure Vessel dimension moddeling, fill and evacuation state determination
V:: 0.0
"""

from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
from scipy import *



def PressureVessel(P, vessel_radius, thickness):
    sigma1 = P*vessel_radius/thickness #calculates hoop stress
    sigma2 = P*vessel_radius/(2*thickness) #calculates longitudinal stress

    print("Hoop Stress :: {:.3e} Pa".format(sigma1))
    print("Longitudinal Stress :: {:.3e} Pa".format(sigma2))
    
    if vessel_radius/thickness < 10:
        print("r/t value too low!")
    #From Matweb.com | Aluminum 6061-T6; 6061-T651
    sigma_yield_T6 = 262000000 #Pa @ 100C 214MPa @149C 103MPa @204C 283MPa @-28C
    #sigma_ultimate_T6 = 290000000 #Pa @100C 234MPa @149 131MPa @204C 324MPa @-28C
    #checks calculated principal stress against selected 
    if sigma1 >= sigma_yield_T6 and sigma2 >= sigma_yield_T6:
        print("Potential Faliure :: Check Hoop & Logitudinal Stresses")
    elif sigma1 >= sigma_yield_T6:
        print("Potential Faliure :: Check Hoop Stress")
    elif sigma2 >= sigma_yield_T6: 
        print("Potential Faliure :: Check Logitudinal Stress")

def calcs():
    T0 = 26.1 #C average temp in June @ NM, Assume N2O @daily mean temp
    vessel_radius = 0.1397 #m | inner radius of vessel
    thickness = 0.00635 #m | vessel wall thickness
    P = 6.895e+6 #Pa | vessel pressure 
    m0 = 20 #kg initial mass
    #Calculates N2O density at given temperature and pressure
    rho_N2O = PropsSI('D', 'T', T0 + 273.15 , 'P', P, 'N2O')
    print("N2O Density {:.3f} kg/m^3".format(rho_N2O))
    
    CheckStress = True
    if CheckStress == True:
        PressureVessel(P, vessel_radius, thickness)

    N2OPhase = CP.PhaseSI('P', P, 'T', T0 +273.15, 'N2O')
    
    U_0 = PropsSI('Q', 'T', T0 + 273.15 , 'P', P, 'N2O')
    print("Initial Internal Energy {} J".format(U_0))

    print(N2OPhase)
    m_dot = .234 #kg/s

    t = 0.04 #s 

    test = [[],[]]

    m_tank = m0
    m_out = 0 #no mass out initially
    while abs(m_out) < m_tank:
        m_out += -m_dot * t
        print(m_out)

    test.append(t)
    test.append(m_tank)

calcs()