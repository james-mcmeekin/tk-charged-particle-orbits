import numpy as np
import matplotlib.pyplot as plt
from borisstep import boris_step, init_half_back

#consts (SI)
q = 1.602176634e-19      #proton charge, C
m = 1.67262192e-27       #proton mass, kg
c = 2.99792458e8         #speed of light, m/s

#tokamak parameters
B0 = 5.0                        # magnetic field strength at the center of the tokamak, T, scalar
R0 = 1.7                        # major radius of tokamak (central axis to to center of plasma cross-section)
q_safety = 2.0                  #safety factor, number of toroidal field lines per poloidal field lines                 
E = np.array([0.0, 0.0, 0.0])

def B_tokamak(pos, B0=B0, R0=R0, q_safety = q_safety):
    x,y,z = pos
    R = np.hypot(x,y)
    #magnetic field in cylindrical coordinates (R, phi, Z)
    B_phi = B0 * R0 / R
    B_R = -B0 * z / (q_safety * R)
    B_z = B0 * (R - R0) / (q_safety * R)

    #convert to Cartesian coordinates (x, y, z)
    B_x = B_R * x / R - B_phi * y / R
    B_y = B_R * y / R + B_phi * x / R

    return np.array([B_x, B_y, B_z])


#launch conditions
v0 = 1.0e6
r0 = 0.1                                         # minor radius of launch point
xi = 0.2                                         # v_|| / v at launch, pitch angle parameter

x = np.array([R0 + r0, 0.0, 0.0], dtype=float)
b_hat = B_tokamak(x, B0=B0, R0=R0, q_safety=q_safety)
b_hat = b_hat / np.linalg.norm(b_hat)

x_hat = np.array([1.0, 0.0, 0.0], dtype=float)


e1 = x_hat - np.dot(x_hat, b_hat) * b_hat       # component of x_hat perpendicular to b_hat
e1 = e1 / np.linalg.norm(e1)
v_init = v0 * (xi * b_hat + np.sqrt(1 - xi**2) * e1)   # initial velocity v(0)


omega_c = B0 * q / m                            # on-axis cyclotron frequency (B0 is a scalar)
T = 2 * np.pi / omega_c
rho = m * v0 * np.sqrt(1 - xi**2) / (q * B0)

steps_per_period = 50
n_periods = 20000                # ~ a few bounce periods for a trapped particle

dt = T / steps_per_period
nsteps = steps_per_period * n_periods

# u(-dt/2), using the field at x(0), so the first boris_step gives u(+dt/2) and x(dt)
u = init_half_back(x, v_init, E, B_tokamak(x, B0=B0, R0=R0, q_safety=q_safety), dt)

traj = np.zeros((nsteps, 3))
u_traj = np.zeros((nsteps, 3))

for i in range(nsteps):
    B = B_tokamak(x, B0=B0, R0=R0, q_safety=q_safety)
    x, u = boris_step(x, u, E, B, dt)
    traj[i], u_traj[i] = x, u


#  diagnostics 
t_gyro = np.arange(1, nsteps + 1) * dt / T

R = np.hypot(traj[:,0], traj[:,1])
z = traj[:,2]
r_minor = np.hypot(R - R0, z)                       # minor radius
theta = np.arctan2(z, R - R0)                       # poloidal angle

# fields along the orbit
B_traj = np.array([B_tokamak(p, B0=B0, R0=R0, q_safety=q_safety) for p in traj])
Bmag = np.linalg.norm(B_traj, axis=1)
v_par = np.sum(u_traj * B_traj, axis=1) / Bmag
u2 = np.sum(u_traj**2, axis=1)
v_perp2 = u2 - v_par**2

# conserved quantities (gamma ~ 1)
KE_err = np.abs(u2 - u2[0]) / u2[0]                                  # energy
mu = m * v_perp2 / (2 * Bmag)                                        # magnetic moment (adiabatic invariant)
mu_err = (mu - mu[0]) / mu[0]
L_z = m * (traj[:,0] * u_traj[:,1] - traj[:,1] * u_traj[:,0])        # m R u_phi
P_phi = L_z + q * B0 * (r_minor**2) / (2 * q_safety)                 # canonical toroidal momentum, psi = B0 r^2 / (2 q_safety)
P_err = (P_phi - P_phi[0]) / np.abs(L_z[0])

bounces = np.count_nonzero(np.diff(np.sign(v_par)) != 0)
xi_crit = np.sqrt(2 * r0 / (R0 + r0))                # trapped if |xi| < xi_crit at the outboard launch point
print(f"{nsteps} steps, dt = T/{T/dt:.0f}, gyroradius = {rho:.3e} m, run = {n_periods * T:.3e} s")
print(f"launch xi = {xi}, trapping threshold ~ {xi_crit:.3f} -> expected {'trapped' if abs(xi) < xi_crit else 'passing'}")
print(f"v_par sign changes (mirror points): {bounces}")
print(f"minor radius range: {r_minor.min():.4f} to {r_minor.max():.4f} m  (launch r0 = {r0})")
print(f"max |dKE/KE0|      = {KE_err.max():.3e}")
print(f"max |dmu/mu0|      = {np.abs(mu_err).max():.3e}")
print(f"max |dP_phi|/|Lz0| = {np.abs(P_err).max():.3e}")

fig, ax = plt.subplots(2, 3, figsize=(17, 9))

# poloidal cross-section
ax[0,0].plot(R, z, lw=0.3)
ax[0,0].plot(R[0], z[0], 'go', label='start')
th = np.linspace(0, 2 * np.pi, 200)
ax[0,0].plot(R0 + r0 * np.cos(th), r0 * np.sin(th), 'k--', lw=0.8, label='r = r0')
ax[0,0].plot(R0, 0, 'k+')
ax[0,0].set_aspect('equal'); ax[0,0].set_title('Poloidal cross-section')
ax[0,0].set_xlabel('R [m]'); ax[0,0].set_ylabel('z [m]'); ax[0,0].legend()
 
# top-down view
ax[0,1].plot(traj[:,0], traj[:,1], lw=0.3)
ax[0,1].set_aspect('equal'); ax[0,1].set_title('Top-down view')
ax[0,1].set_xlabel('x [m]'); ax[0,1].set_ylabel('y [m]')

# v_par / v showing mirroring
ax[0,2].plot(t_gyro, v_par / np.sqrt(u2), lw=0.5)
ax[0,2].axhline(0, color='k', lw=0.5)
ax[0,2].set_title('Pitch v_par / v (sign change = trapped)')
ax[0,2].set_xlabel('gyroperiods')

ax[1,0].semilogy(t_gyro, np.maximum(KE_err, 1e-17))
ax[1,0].set_title('Energy error |dKE/KE0|'); ax[1,0].set_xlabel('gyroperiods')

ax[1,1].plot(t_gyro, mu_err, lw=0.5)
ax[1,1].set_title('Magnetic moment error (mu - mu0)/mu0'); ax[1,1].set_xlabel('gyroperiods')

ax[1,2].plot(t_gyro, P_err, lw=0.5)
ax[1,2].set_title('Canonical P_phi error / |Lz0|'); ax[1,2].set_xlabel('gyroperiods')

plt.tight_layout()
plt.show()
