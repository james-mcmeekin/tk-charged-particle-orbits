import numpy as np
import matplotlib.pyplot as plt

#consts (SI)
q = 1.602176634e-19      #proton charge, C
m = 1.67262192e-27       #proton mass, kg
c = 2.99792458e8         #speed of light, m/s

# remade boris step function because import statement will raise TypeError
def boris_step(x, u, E, B, dt):
    h = q * dt / (2 * m)
    u_minus = u + h * E
    gamma = np.sqrt(1 + np.dot(u_minus, u_minus) / c**2)
    t = h * B / gamma
    s = 2 * t / (1 + np.dot(t, t))
    uprime = u_minus + np.cross(u_minus, t)
    u_plus = u_minus + np.cross(uprime, s)
    unew = u_plus + h * E
    gamma_new = np.sqrt(1 + np.dot(unew, unew) / c**2)
    xnew = x + unew * dt / gamma_new
    return xnew, unew


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

omega_c = q * np.linalg.norm(B) / m
T = 2 * np.pi / omega_c
rho = m * v_perp / (q * np.linalg.norm(B))

dt = T / 50
n_gyro = 1000
nsteps = int(n_gyro * T / dt)

x0 = np.zeros(3)
u0 = np.array([0.0, v_perp, 0.0])

# Boris: u lives at half-integer steps (leapfrog), so push u back by dt/2 first
xb, ub = x0.copy(), u0.copy()
_, ub = boris_step(xb, ub, E, B, -dt / 2)
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

print(f"dt = T/{T/dt:.0f}, {n_gyro} gyroperiods, {nsteps} steps")
print(f"final |dE/E|:  Boris = {dE_boris[-1]:.3e}   RK4 = {dE_rk4[-1]:.3e}")
print(f"final radius error (r-rho)/rho:  Boris = {(r_boris[-1]-rho)/rho:.3e}   RK4 = {(r_rk4[-1]-rho)/rho:.3e}")

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))


# plots position for first 5 gyroperiods
n_show = int(5 * T / dt)
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
