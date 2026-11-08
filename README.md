This project is a simulated boris pusher for charged particles, with results validated with expected results. The code is written in Python and uses NumPy for numerical computations. It simulates the motion of charged particles under the influence of electric and magnetic fields using the Boris algorithm, which is a widely used method for integrating the equations of motion for charged particles.
The results of the simulation are visualized using Matplotlib, and are validated against expected results to ensure accuracy, as well as comparison with the Runge-Kutta method. 


In nuclear fusion, the boris pusher algorithm is used to simulate the charged particles within plasma, and is incredibly effective as it conserves |u|, the magnitude of the proper velocity, through rotation exactly, making it suitable for long-term simulations.

The algorithm works by splitting the motion of the particles into two parts: the acceleration from the magnetic field, and the acceleration from the electric field. 
It first updates the velocity of the particles with a "half-kick" from the electric field, then applies a rotation due to the magnetic field, and finally updates the velocity again with another "half-kick" from the electric field. This method ensures that the particles' motion is accurately represented over time. 

The equation that defines the motion of a charged particle in an electromagnetic field is given by the Lorentz force law:

$$
\begin{aligned}
   \mathbf{F} = q(\mathbf{E} + \mathbf{v} \times \mathbf{B}) 
\end{aligned}$$

Where:
-  $\mathbf{F}$ is the force acting on the particle, measured in newtons (N),
- $q$ is the charge of the particle, measured in coulombs (C),
- $\mathbf{E}$ is the electric field, measured in volts per meter (V/m),
- $\mathbf{v}$ is the velocity of the particle, measured in meters per second (m/s),
- $\mathbf{B}$ is the magnetic field, measured in Tesla (T).


In this algorithm, we define the velocity of the particle as $\mathbf{u} = \gamma \mathbf{v}$, otherwise known as the "proper velocity" where $\gamma$ is the relativistic gamma factor given by:

$$
\begin{aligned}
    \gamma &= \frac{1}{\sqrt{1-\frac{\mathbf{v}^2}{c^2}}} \\
           &= \sqrt{1+\frac{\mathbf{u}^2}{c^2}}
\end{aligned}
$$

Therefore we can write a relativistic set of equations of motion:

$$\begin{aligned}
    \frac{d\mathbf{u}}{dt} = \frac{q}{m}(\mathbf{E} + \frac{\mathbf{u}}{\gamma} \times \mathbf{B}) \\
    \frac{d\mathbf{x}}{dt} = \frac{\mathbf{u}}{\gamma}
\end{aligned}$$


For each step over time $\Delta t$, the algorithm takes in the initial proper velocity $\mathbf{u}$ and performs the first half-kick from the electric field:

$$
\begin{aligned}
    \mathbf{u}^{-} = \mathbf{u} + \frac{q\Delta t}{2m}\mathbf{E}
\end{aligned}$$

It then calculates the intermediate vector $\mathbf{u}'$ that is used to calculate the full rotation:

$$
\begin{aligned}
    \mathbf{t} = \frac{q\Delta t}{2m\gamma^-}\mathbf{B} \\
    \mathbf{u}' = \mathbf{u}^{-} + \mathbf{u}^{-} \times \mathbf{t}
\end{aligned}
$$

Where $\gamma^{-} = \sqrt{1+\frac{(\mathbf{u}^{-})^2}{c^2}}$.

Next, Let us define $\mathbf{u}^{+}$ as the final proper velocity after the rotation.
To calculate the second half of the rotation, the second rotation factor $\mathbf{s}$ found by the rotation condition:

$$
\begin{aligned}
 \mathbf{u}^{+} - \mathbf{u}^{-} = (\mathbf{u}^{+} + \mathbf{u}^{-}) \times \mathbf{t}
\end{aligned}$$

We know that the second half of the rotation is given by:

$$
\begin{aligned}
    \mathbf{u}^{+} = \mathbf{u}^{-} + \mathbf{u}' \times \mathbf{s}
\end{aligned}$$

To solve for $\mathbf{s}$, we first assume $\mathbf{s} = \lambda \mathbf{t}$, where $\lambda$ is a scalar.
Substituting $u^{+}$ into the rotation condition, we have:

$$
\begin{aligned}
    \mathbf{u}^{-} + \mathbf{u}' \times \lambda \mathbf{t} - \mathbf{u}^{-} = (\mathbf{u}^{-} + \mathbf{u}' \times \lambda \mathbf{t} + \mathbf{u}^{-}) \times \mathbf{t} \\
    \mathbf{u}' \times \lambda \mathbf{t} = (2\mathbf{u}^{-} + \mathbf{u}' \times \lambda \mathbf{t}) \times \mathbf{t}
\end{aligned}$$

First solving for $\mathbf{u}' \times \lambda\mathbf{t}$:

$$
\begin{aligned}
    \mathbf{u}' \times \mathbf{s} &= \mathbf{u}^{-} \times \lambda\mathbf{t} + (\mathbf{u}^{-} \times \mathbf{t}) \times \lambda\mathbf{t} \\
    &= \lambda(\mathbf{u}^{-} \times \mathbf{t} + \mathbf{u}^{-} \cdot \mathbf{t})\mathbf{t} - (\mathbf{t}^{2}\mathbf{u}^{-}) \\
\end{aligned}$$

Let us call 
$\mathbf{A} = \mathbf{u}^{-} \times \mathbf{t} + (\mathbf{u}^{-} \cdot \mathbf{t})\mathbf{t} - \mathbf{t}^{2}\mathbf{u}^{-}$
, then we have:

$$
\begin{aligned}
    \mathbf{u}' \times \lambda\mathbf{t} = \mathbf{A}\lambda
\end{aligned}$$
    

Solving the right hand side of the rotation condition:

$$
\begin{aligned}
    (2\mathbf{u}^{-}+\mathbf{A}\lambda) \times \mathbf{t} &= 2\mathbf{u}^{-} \times \mathbf{t} + \mathbf{A}\lambda \times \mathbf{t} \\
    &= 2\mathbf{u}^{-} \times \mathbf{t} + \lambda((\mathbf{u}^{-}\times \mathbf{t}) \times \mathbf{t} + (\mathbf{u}^{-} \cdot \mathbf{t})\mathbf{t} \times \mathbf{t} - (\mathbf{t}^{2}\mathbf{u}^{-}\times \mathbf{t})) \\
    &= 2\mathbf{u}^{-} \times \mathbf{t} + \lambda((\mathbf{u}^{-}\times \mathbf{t}) \times \mathbf{t} - (\mathbf{t}^{2}\mathbf{u}^{-}\times \mathbf{t})) \\
\end{aligned}$$

Therefore we have:

$$
\begin{aligned}
    \mathbf{A}\lambda = 2\mathbf{u}^{-} \times \mathbf{t} + \lambda((\mathbf{u}^{-}\times \mathbf{t}) \times \mathbf{t} - (\mathbf{t}^{2}(\mathbf{u}^{-}\times \mathbf{t}))) \\
    \lambda(\mathbf{u}^{-} \times \mathbf{t} + (\mathbf{u}^{-} \cdot \mathbf{t})\mathbf{t} - \mathbf{t}^{2}\mathbf{u}^{-}) = 2\mathbf{u}^{-} \times \mathbf{t} + \lambda((\mathbf{u}^{-}\times \mathbf{t}) \times \mathbf{t} - \mathbf{t}^{2}(\mathbf{u}^{-}\times \mathbf{t})) \\
    \lambda(\mathbf{u}^{-} \times \mathbf{t} + (\mathbf{u}^{-} \cdot \mathbf{t})\mathbf{t} - \mathbf{t}^{2}\mathbf{u}^{-}) = 2\mathbf{u}^{-} \times \mathbf{t} + \lambda((\mathbf{u}^{-}\cdot \mathbf{t})\mathbf{t} - \mathbf{t}^{2}\mathbf{u}^{-} - \mathbf{t^2}(\mathbf{u}^{-}\times \mathbf{t}))
\end{aligned}$$
    
Canceling equal terms:

$$
\begin{aligned}
    \lambda(\mathbf{u}^{-} \times \mathbf{t}) = 2(\mathbf{u}^{-} \times \mathbf{t}) - \lambda(t^2\mathbf{u}^{-} \times \mathbf{t}) \\
    \implies \lambda(t^2+1)=2 \\
    \implies \lambda = \frac{2}{1+t^2}
\end{aligned}$$

Therefore

$$
\begin{aligned}
    \mathbf{s} = \frac{2\mathbf{t}}{1+t^2}
\end{aligned}$$
    
Therefore, returning to the boris step function, we calculate the second half of the rotation as:

$$
\begin{aligned}
    \mathbf{u}^{+} = \mathbf{u}^{-} + \mathbf{u}' \times \mathbf{s}
\end{aligned}
$$
    
And finally, it performs the second half-kick from the electric field to update the proper velocity:

$$
\begin{aligned}
    \mathbf{u}_{new} = \mathbf{u}^{+} + \frac{q\Delta t}{2m}\mathbf{E}
\end{aligned}$$

To put each step in order, the algorithm is as follows:
1. Calculate the first half-kick from the electric field:
   $\mathbf{u}^{-} = \mathbf{u} + \frac{q\Delta t}{2m}\mathbf{E}$
2. Calculate the rotation vectors:
   $\mathbf{t} = \frac{q\Delta t}{2m\gamma^-}\mathbf{B}$
   ,
   $\mathbf{s} = \frac{2\mathbf{t}}{1+\mathbf{t}^2}$
3. Calculate the intermediate rotation vector:
   $\mathbf{u}' = \mathbf{u}^{-} + \mathbf{u}^{-} \times \mathbf{t}$
4. Calculate the rotation:
   $\mathbf{u}^{+} = \mathbf{u}^{-} + \mathbf{u}' \times \mathbf{s}$
5. Calculate the second half-kick from the electric field:
   $\mathbf{u}_{new} = \mathbf{u}^{+} + \frac{q \Delta t}{2m}\mathbf{E}$

Which we can update the position of the particle as:

$$
\begin{aligned}
    \mathbf{x}_{new} = \mathbf{x} + \frac{\mathbf{u}_{new}}{\gamma_{new}}\Delta t
\end{aligned}
$$

We can also see how $\mathbf{|u|}$ is conserved through the rotation algorithm:

$$
\begin{aligned}
    \mathbf{u}^{+}-\mathbf{u}^{-} = (\mathbf{u}^{+} + \mathbf{u}^{-}) \times \mathbf{t}
    (\mathbf{u}^{+}-\mathbf{u}^{-}) \cdot (\mathbf{u}^{+} + \mathbf{u}^{-}) = ((\mathbf{u}^{+} + \mathbf{u}^{-}) \times \mathbf{t}) \cdot (\mathbf{u}^{+} + \mathbf{u}^{-}) \\
\end{aligned}
$$

On the right hand side, the cross product term is by definition perpendicular to $(\mathbf{u}^{+} + \mathbf{u}^{-})$, and therefore:

$$
\begin{aligned}
    (\mathbf{u}^{+}-\mathbf{u}^{-}) \cdot (\mathbf{u}^{+} + \mathbf{u}^{-}) = 0 \\
    \implies |\mathbf{u}^{+}|^2 - |\mathbf{u}^{-}|^2 = 0 \\
    \implies |\mathbf{u}^{+}| = |\mathbf{u}^{-}|
\end{aligned}
$$

So magnitude is clearly preserved through rotation. Magnitude will therefore only change in the presence of an electric field (i.e $\mathbf{E}\neq 0$), which is why we calculate a new gamma factor after the second half-kick to update the position correctly.

One caveat to the boris algorithm, is that the rotation angle is calculated as:

$$
\begin{aligned}
    \theta = 2 \arctan(|\mathbf{t}|)
\end{aligned}
$$

Rather than $\frac{\omega_c \cdot \Delta t}{\gamma}$, where $\omega_c = \frac{q|\mathbf{B}|}{m}$, denoted as the Cyclotron frequency, and that error is the source of the bounded error we see when comparing to RK4.

Implementation Notes to borisstep.py

For accurate position updating, we have to adjust $\mathbf{u}$ by a half-step backwards in time before starting borisstep. This is because we want the velocity at the midpoint between $\Delta t_n$ and $\Delta t_{n+1}$, meaning to adjust $\mathbf{x}_{0} \to \mathbf{x}_{1}$, we want $\mathbf{u}_{\frac{1}{2}}$.
That is, to update the position correctly, we actually compute:

$$
\begin{aligned}
    \mathbf{x}_{n+1} = \mathbf{x}_{n}+ \frac{\mathbf{u}_{n+\frac{1}{2}}}{\gamma_{n+\frac{1}{2}}}\Delta t
\end{aligned}
$$

So to initially get $\mathbf{u}_{-\frac{1}{2}}$, we call the `borisstep` function with a negative half-step ($-\frac{\Delta t}{2}$) in time via init_half_back(). If this weren't implemented, then we would have a startup error of $O(\delta t)$, which would take precedent over the algorithms actual error of $O(\delta t^2)$.

The plot generated by `borisstep.py` shows the particle trajectory in the x-y plane on the left and the energy error relative to the initial kinetic energy on the right. For constant a constant magnetic field, the energy error should remain small, which is verified by the plot being on the order of ~$10^{-17}$-$10^{-15}$.

Boris Step vs. RK4

 `compare_rk4.py` compares the Boris step algorithm with the RK4 method for a charged particle in a constant magnetic field. The script generates plots of the particle trajectories, energy errors and gyroradius errors for both methods, allowing for a direct comparison of accuracy and stability.

 For a constant magnetic field, we define the accepted gyroradius as:

$$
\rho = \frac{\gamma v_\perp}{\omega_c}
$$

where $\omega_c = \frac{q|\mathbf{B}|}{m}$ is the cyclotron frequency. This simply stems from the equation between angular velocity and linear velocity for circular motion, $v_\perp = \omega_c \rho$.

In order to calculate the error in the gyroradius, we need to compare the position of the particle relative to the guiding center of its orbit. We know the particle will orbit around a circle of radius rho, and we can find that center through the initial conditions at $t=0$:

$$
\begin{aligned}
    \mathbf{v}_\perp = v_\perp \cdot \hat{y} \\
    \mathbf{B} = B \cdot \hat{z}
\end{aligned}
$$

Calculating the force, which points towards the center of the circular orbit, we get:

$$
\begin{aligned}
\mathbf{F} &= q \mathbf{v}_\perp \times \mathbf{B} \\
&= q \mathbf{v}_\perp \mathbf{B} (\hat{y} \times \hat{z}) \\
&= q \mathbf{v}_\perp \mathbf{B} \hat{x}
\end{aligned}
$$

Therefore, the guiding center of the particle's orbit is located at $(\rho, 0, 0)$. So in order to calculate the error in the gyroradius, we compare the particle's position relative to this guiding center;

$$
\begin{aligned}
    \mathbf{r} &= \sqrt{(x-\rho)^2 + y^2 + z^2} 
\end{aligned}
$$

Which since the particle in this case does not move in the $z$ direction, simplifies to:

$$
\begin{aligned}
    \mathbf{r} &= \sqrt{(x-\rho)^2 + y^2} \\
\end{aligned}
$$


 We see from the plot that the Boris Algorithm has a constant range of energy error independent of time, while RK4 starts off with a smaller error but gradually increases over time, which is why the Boris algorithm is preferred for long-term simulations of charged particle motion in magnetic fields. Plotting over longer gyroperiods (e.g 50000) will make the difference apparent
