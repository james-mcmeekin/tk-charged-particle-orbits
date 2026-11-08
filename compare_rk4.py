import numpy as np
import matplotlib.pyplot as plt
from borisstep import boris_step, init_half_back

#consts (SI)
q = 1.602176634e-19      #proton charge, C
m = 1.67262192e-27       #proton mass, kg
c = 2.99792458e8         #speed of light, m/s

#Right hand side of Runge-Kutta method, dy/dt = f(y, t); ignore t variable because it is not used
def rhs(y, E, B):
    # y = [x, u]; dx/dt = u/gamma, du/dt = q/m (E + v x B)
    u = y[3:]
    gamma = np.sqrt(1 + np.dot(u, u) / c**2)
    v = u / gamma
    return np.concatenate((v, q / m * (E + np.cross(v, B))))

#runge-kutta 4th order stepper
def rk4_step(y, E, B, dt):
    k1 = rhs(y, E, B)
    k2 = rhs(y + 0.5 * dt * k1, E, B)
    k3 = rhs(y + 0.5 * dt * k2, E, B)
    k4 = rhs(y + dt * k3, E, B)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

#relativistic kinetic energy function KE = (1-gamma)mc^2
def kinetic_energy(u):
    gamma = np.sqrt(1 + np.sum(u * u, axis=-1) / c**2)
    return (gamma - 1) * m * c**2

# setup ics (same as borisstep.py)
B = np.array([0.0, 0.0, 5.0])
E = np.array([0.0, 0.0, 0.0])
v_perp = 1.0e6
gamma = 1 / np.sqrt(1 - v_perp**2 / c**2)

omega_c = q * np.linalg.norm(B) / m
T = 2 * np.pi / omega_c
rho = gamma * v_perp / omega_c

steps_per_period = 50
n_periods = 1000

dt = T / steps_per_period
nsteps = steps_per_period * n_periods

x0 = np.zeros(3)
v0 = np.array([0.0, v_perp, 0.0])
u0 = v0 / np.sqrt(1 - np.dot(v0, v0) / c**2)   # u = gamma * v, shared by Boris and RK4

# Boris: u lives at half-integer steps (leapfrog), so push u back by dt/2 first
xb = x0.copy()
ub = init_half_back(xb, v0, E, B, dt)   # u(-dt/2) from v(0)
boris_x = np.zeros((nsteps, 3))
boris_u = np.zeros((nsteps, 3))
for i in range(nsteps):
    xb, ub = boris_step(xb, ub, E, B, dt)
    boris_x[i], boris_u[i] = xb, ub

# RK4, same dt
# y is of form [x, y, z, u_x, u_y, u_z]
y = np.concatenate((x0, u0))
rk_x = np.zeros((nsteps, 3))
rk_u = np.zeros((nsteps, 3))
for i in range(nsteps):
    y = rk4_step(y, E, B, dt)
    rk_x[i], rk_u[i] = y[:3], y[3:]

t_gyro = (np.arange(1, nsteps + 1) * dt) / T
KE0 = kinetic_energy(u0)

#calculate relative kinetic energy error for boris and rk4 methods
dE_boris = np.abs(kinetic_energy(boris_u) - KE0) / KE0
dE_rk4 = np.abs(kinetic_energy(rk_u) - KE0) / KE0


r_boris = np.hypot(boris_x[:, 0] - rho, boris_x[:, 1])
r_rk4 = np.hypot(rk_x[:, 0] - rho, rk_x[:, 1])

print(f"dt = T/{T/dt:.0f}, {n_periods} gyroperiods, {nsteps} steps")
print(f"final |dE/E|:  Boris = {dE_boris[-1]:.3e}   RK4 = {dE_rk4[-1]:.3e}")
print(f"final radius error (r-rho)/rho:  Boris = {(r_boris[-1]-rho)/rho:.3e}   RK4 = {(r_rk4[-1]-rho)/rho:.3e}")

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))


# plots position for first 5 gyroperiods
n_show = round(5 * steps_per_period)
ax[0].plot(boris_x[:n_show, 0], boris_x[:n_show, 1], label='Boris')
ax[0].plot(rk_x[:n_show, 0], rk_x[:n_show, 1], '--', label='RK4')
ax[0].set_aspect('equal'); ax[0].set_title('First 5 gyroperiods')
ax[0].set_xlabel('x [m]'); ax[0].set_ylabel('y [m]'); ax[0].legend()

#plots energy error over time (gyroperiods); boris should have energy error of 0.
ax[1].semilogy(t_gyro, dE_boris, label='Boris')
ax[1].semilogy(t_gyro, dE_rk4, label='RK4')
ax[1].set_title('Energy error |ΔKE/KE₀|')
ax[1].set_xlabel('gyroperiods'); ax[1].legend()


#plots gyroradius error over time (gyroperiods). Boris should have constant range of error, 
# while RK4 should increasingly drift over time.
ax[2].plot(t_gyro, (r_boris - rho) / rho, label='Boris')
ax[2].plot(t_gyro, (r_rk4 - rho) / rho, label='RK4')
ax[2].set_title('Gyroradius error (r-ρ)/ρ')
ax[2].set_xlabel('gyroperiods'); ax[2].legend()

plt.tight_layout()
plt.show()
