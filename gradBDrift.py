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

vperp = 1e6
rho = m * vperp / (q * np.linalg.norm(B0)) # ~ 2 mm
L = 100 * rho                              # gradient scale length

omega_c = np.linalg.norm(B) * q / m
T = 2 * np.pi / omega_c
dt = T / 50
n_gyro = 1000
