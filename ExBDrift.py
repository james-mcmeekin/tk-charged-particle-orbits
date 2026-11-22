import numpy as np  
import matplotlib.pyplot as plt
from borisstep import boris_step, init_half_back

#consts (SI)
q = 1.602176634e-19      #proton charge, C
m = 1.67262192e-27     #proton mass, kg

#ics:
B = np.array([0,0,5], dtype=float)     #Magnetic Field (T)
E = np.array([1e4,0,0], dtype=float)   #Electric Field (V/m)
v_perp = 1.0e6

omega_c = q * np.linalg.norm(B) / m          # nonrelativistic cyclotron frequency
T = 2 * np.pi / omega_c                        # cyclotron gyroperiod
rho = m * v_perp / (q * np.linalg.norm(B))   # analytic gyroradius

n_periods = 20
steps_per_period = 50

dt = T / steps_per_period
nsteps = steps_per_period * n_periods

x = np.array([0.0, 0.0, 0.0], dtype = float)       # initial position
v0 = np.array([0.0, v_perp, 0.0], dtype = float)   # initial velocity v(0)
u = init_half_back(x, v0, E, B, dt)                # u(-dt/2), so the first boris_step gives u(+dt/2) and x(dt)

traj = np.zeros((nsteps + 1, 3))                   # traj[0] is the starting point
traj[0] = x
for i in range(nsteps):
    x,u = boris_step(x, u, E, B, dt)
    traj[i + 1] = x

t = dt * np.arange(nsteps + 1)

# theoretical charge drift: v_d = E x B / B^2 = E/B * (xhat cross zhat) = E/B yhat = 2000 m/s in y direction
v_drift_theory = np.cross(E, B) / np.dot(B, B)

# measured: fit pos(t) = x0 + v_d*t + gyration. Boris rotates by 2*arctan(omega_c*dt/2) per step,
# so the gyration uses this numerical frequency (a plain endpoint difference (e.g u from boris_step) would leak phase error into the drift)
omega_num = 2 * np.arctan(omega_c * dt / 2) / dt
A = np.column_stack([np.ones_like(t), t, np.cos(omega_num * t), np.sin(omega_num * t)])
coef, *_ = np.linalg.lstsq(A, traj, rcond=None)    # one column per x,y,z
v_drift_sim = coef[1]
gc0 = coef[0]                                      # guiding center

print(f"theory    v_d = {v_drift_theory} m/s")
print(f"simulated v_d = {v_drift_sim} m/s")
print(f"relative error = {np.linalg.norm(v_drift_sim - v_drift_theory) / np.linalg.norm(v_drift_theory):.3e}")

# guiding center path = trajectory with the theoretical drift removed should be a closed circle
gc_frame = traj - np.outer(t, v_drift_theory)
r = np.hypot(gc_frame[:,0] - gc_frame[:,0].mean(), gc_frame[:,1] - gc_frame[:,1].mean())
print(f"gyroradius in drift frame: analytic = {rho:.3e} m, simulated mean = {r.mean():.3e} m")

#plot orbit of particle
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
ax[0].plot(traj[:,0], traj[:,1])
ax[0].set_xlabel('x [m]'); ax[0].set_ylabel('y [m]')
ax[0].set_title('Orbit in lab frame')
ax[0].set_aspect('equal', adjustable='datalim')

#plots of the particle's y-coordinate and the guiding center's y-coordinate over time
ax[1].plot(t / T, traj[:,1], label='simulated y(t)')
ax[1].plot(t / T, gc0[1] + v_drift_theory[1] * t, '--', label='theory guiding center')
ax[1].set_xlabel('t / T_c'); ax[1].set_ylabel('y [m]')
ax[1].set_title('Drift along E x B (-y)')
ax[1].legend()
plt.tight_layout()
plt.show()
