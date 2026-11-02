import numpy as np
import matplotlib.pyplot as plt
from borisstep import boris_step

#consts (SI)
q = 1.602176634e-19      #proton charge, C
m = 1.67262192e-27       #proton mass, kg
c = 2.99792458e8         # speed of light, m/s


#ics:
B = np.array([0.0, 0.0, 5.0], dtype=float) # Tesla                                  # gradient scale length
E = np.array([0.0, 0.0, 0.0], dtype=float)  # V/m

def B_field(pos, B=B, L=0.2):
    y = pos[1]
    return np.array([B[0], B[1], B[2] * (1 + y / L)], dtype=float)

v_perp = 1e6


omega_c = np.linalg.norm(B) * q / m
T = 2 * np.pi / omega_c
rho = m * v_perp / (q * np.linalg.norm(B)) 
L = 100 * rho                               #gradient scale length

steps_per_period = 50
n_periods = 1000

dt = T / steps_per_period
nsteps = steps_per_period * n_periods
