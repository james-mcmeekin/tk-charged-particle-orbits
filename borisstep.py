import numpy as np
import matplotlib.pyplot as plt


#consts (SI)
q = 1.602176634e-19      #proton charge, C
m = 1.67262192e-27       #proton mass, kg
c = 2.99792458e8         #speed of light, m/s

def boris_step(x, u, E, B,dt):  
    h = q * dt / (2 * m)                                  # half electric kick
    u_minus = u + h * E                                   # first intial half electric kick
    gamma = np.sqrt(1 + np.dot(u_minus, u_minus) / c**2)  # relativistic gamma factor
    t = h * B / gamma                                     
    s = 2 * t / (1 + np.dot(t, t)) 
    uprime = u_minus + np.cross(u_minus, t)               # rotation intermediate vector
    u_plus = u_minus + np.cross(uprime, s)                # rotation
    unew = u_plus + h * E                                 # second half electric kick
    gamma_new = np.sqrt(1 + np.dot(unew, unew) / c**2)      
    xnew = x + unew * dt / gamma_new                      # new pos
    return xnew, unew

def init_half_back(x0, v0, E, B, dt):
    #given v(0) (ordinary velocity), return u(-dt/2) so that boris_step(x0, u(-dt/2), ...) yields u(+dt/2), x(dt).
    gamma0 = 1 / np.sqrt(1 - np.dot(v0, v0) / c**2)  # gamma from ordinary velocity v (not u)
    u0 = v0 * gamma0                                 # u = gamma * v
    _, uback = boris_step(x0, u0, E, B, -dt / 2)     # Boris is time-reversible: negative dt/2 steps u back by dt/2; position output discarded
    return uback



#ensures the boris step experiment is ran only when the file is ran; boris_step used in other files
if __name__ == "__main__":
    # setup ics
    B = np.array([0.0, 0.0, 5.0], dtype = float)     # Tesla
    E = np.array([0.0, 0.0, 0.0], dtype = float)     # V/m
    v_perp = 1.0e6                                   # m/s
    gamma = 1 / np.sqrt(1 - v_perp**2 / c**2)        # initial relativistic gamma factor

    omega_c = q * np.linalg.norm(B) / m                   # cyclotron frequency; how many times a particle "orbits" magnetic field line.
    T = 2 * np.pi / omega_c                               # cyclotron gyroperiod; period of one full orbit around the magnetic field line.
    rho = gamma * v_perp / omega_c                        # analytic (accepted) gyroradius; rho = gamma*m*v_perp / (|q|*B)

    n_periods = 5
    steps_per_period = 50

    dt = T / steps_per_period
    nsteps = steps_per_period * n_periods

    x = np.array([0.0, 0.0, 0.0], dtype = float)               # initial position
    v0 = np.array([0.0, v_perp, 0.0], dtype = float)           # intial velocity
    u = init_half_back(x, v0, E, B, dt)                        # initial true velocity, adjusted for boris step

    traj = np.zeros((nsteps, 3))
    u_traj = np.zeros((nsteps, 3))
    for i in range(nsteps):
        x,u = boris_step(x, u, E, B, dt)
        traj[i] = x
        u_traj[i] = u

    KE = np.sum(u_traj**2, axis=1)
    dE = np.abs(KE - KE[0]) / KE[0]

    radius = np.hypot(traj[:,0] - rho, traj[:,1])
    print(f"analytic rho = {rho:.3e} m")
    print(f"simulated radius: min = {np.min(radius):.3e} m, max = {np.max(radius):.3e} m, error = {np.max(radius) - np.min(radius):.3e} m")

    fig, ax = plt.subplots(1, 2, figsize=(16, 4.5))

    ax[0].plot(traj[:,0], traj[:,1])
    ax[0].set_aspect('equal')
    ax[0].set_title('Particle trajectory in x-y plane')
    ax[0].set_xlabel('x [m]')
    ax[0].set_ylabel('y [m]')

    t_gyro = np.arange(1, nsteps + 1) * dt / T
    ax[1].semilogy(t_gyro, np.maximum(dE, 1e-17))
    ax[1].set_title('Energy error |ΔKE/KE₀|')
    ax[1].set_xlabel('gyroperiods')

    plt.tight_layout()
    plt.show()
