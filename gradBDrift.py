import numpy as np
import matplotlib.pyplot as plt
from borisstep import boris_step, init_half_back

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
v0 = np.array([v_perp, 0.0, 0.0], dtype=float)   # initial velocity v(0)
u = init_half_back(x, v0, E, B_field(x, B0, L), dt)   # u(-dt/2), using the field at x(0)

traj = np.zeros((nsteps, 3))
u_traj = np.zeros((nsteps, 3))
for i in range(nsteps):
    B = B_field(x, B0, L)   # Bz = B0[2] * (1 + y / L)
    x, u = boris_step(x, u, E, B, dt)
    traj[i], u_traj[i] = x, u

x_avg = traj[:,0].reshape(n_periods, steps_per_period).mean(axis=1)
t_avg = np.arange(n_periods) * T
slope = np.polyfit(t_avg, x_avg, 1)[0]

v_pred = -m * v_perp**2 / (2 * q * B0[2] * L)
print(f"simulated drift: {slope:.1f} m/s, predicted: {v_pred:.1f} m/s")

# energy check: B does no work, so |u| (and KE) should be conserved
KE = np.sum(u_traj**2, axis=1)
dE = np.abs(KE - KE[0]) / KE[0]

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))

# trajectory for first 30 gyroperiods
n_show = 30 * steps_per_period
ax[0].plot(traj[:n_show, 0], traj[:n_show, 1])
ax[0].set_aspect('equal'); ax[0].set_title('First 30 gyroperiods')
ax[0].set_xlabel('x [m]'); ax[0].set_ylabel('y [m]')

# period-averaged x vs time with predicted drift
ax[1].plot(t_avg / T, x_avg, '.', label='simulated (period avg)')
ax[1].plot(t_avg / T, x_avg[0] + v_pred * t_avg, '--', label=f'predicted ({v_pred:.0f} m/s)')
ax[1].set_title('Guiding-centre x vs time')
ax[1].set_xlabel('gyroperiods'); ax[1].set_ylabel('<x> [m]'); ax[1].legend()

t_gyro = np.arange(1, nsteps + 1) * dt / T
ax[2].semilogy(t_gyro, np.maximum(dE, 1e-17))
ax[2].set_title('Energy error |ΔKE/KE₀|')
ax[2].set_xlabel('gyroperiods')

plt.tight_layout()
plt.show()
