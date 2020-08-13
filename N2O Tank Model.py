"""
OOrtiz
01AUG20
V:: 0.0
"""

from sympy import *
import numpy as np
from sympy import symbols, Eq, solve
from sympy.solvers import solve
from sympy import Symbol
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def model(z,t):
    k = 0.3
    dydt = -k * y
    dxdt = 1
    dzdt = [dydt, dxdt]
    return dzdt

y0 = 5
z0 = [0,0]
t = np.linspace(0,20)

y = odeint(model, z0, t)

plt.plot( t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()



#dT_wall/dt = x(1)
#dT_wall/dt = (Qdot_wall_in - Qdot_wall_out + Q_dot_wall_cond + m_dot_wall_in * c_wall * (T_wall_in - T_wall)) / (m_wall * c_wall)
#dQ_wall_cond/dt = x(2) 
#dQ_wall_cond = (k_wall * (Delta T)_wall * pi * (r_0**2 - r_i**2))/L_wal_cond
#   L_wall_cond = distance between liquid and vapor sections of the wall



# T_wall_liq = x(1)
# dt_wall_vap = x(2)

###used in x(1)
# Q_wall_in_liq = x(3)
# Q_wall_out_liq = x(4)
# Q_wall_cond_liq = x(5)
# m_wall_in_liq = x(6)

### used in x(2)
# Q_wall_in_vap = x(7)
# Q_wall_out_vap = x(8)
# Q_wall_cond_vap = x(9)
# m_wall_in_vap = x(10)

# m_vap = x(11)
# m_liq = x(12)

### used in x(11)
#m_evap = x(13) ###Used in x(12) too
#m_condensate = x(14)

### used in x(12)
#m_outlet = x(13)

# U_vap = x(14)
# U_liq = x(15)

### used in x(14)
# V_vap = x(16)
# Q_in_vap = x(17)

### used in x(15)
# V_liq = x(18)
# Q_in_liq = x(19)

# T_liq = x(20)
# T_vap = x(21)