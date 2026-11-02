import numpy as np
from borisstep import boris_step

#consts (SI)
q = 1.602176634e-19      #proton charge, C
m = 1.67262192e-27       #proton mass, kg


#ics:
B0 = np.array([0.0, 0.0, 5.0], dtype=float) # Tesla
E = np.array([0.0, 0.0, 0.0], dtype=float)  # V/m

def B_field(pos, B, L):
    y = pos[1]
    return np.array([B[0], B[1], B[2] * (1 + y / L)], dtype=float)

v_perp = 1e6


omega_c = np.linalg.norm(B0) * q / m
T = 2 * np.pi / omega_c
rho = m * v_perp / (q * np.linalg.norm(B0))
L = 100 * rho                               # gradient scale length

steps_per_period = 50
n_periods = 1000

dt = T / steps_per_period
nsteps = steps_per_period * n_periods

# position and velocity such that the particles guiding center is at origin
x = np.array([0.0, rho, 0.0], dtype=float)       # initial position
u = np.array([v_perp, 0.0, 0.0], dtype=float)    # initial velocity, allowing gamma ~ 1

traj = np.zeros((nsteps, 3))
for i in range(nsteps):
    B = B_field(x, B0, L)   # Bz = B0[2] * (1 + y / L)
    x, u = boris_step(x, u, E, B, dt)
    traj[i] = x

x_avg = traj[:,0].reshape(n_periods, steps_per_period).mean(axis=1)
t_avg = np.arange(n_periods) * T
slope = np.polyfit(t_avg, x_avg, 1)[0]

v_pred = -m * v_perp**2 / (2 * q * B0[2] * L)
print(f"simulated drift: {slope:.1f} m/s, predicted: {v_pred:.1f} m/s")
