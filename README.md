This project is a simulated boris pusher for charged particles, with results validated with expected results. The code is written in Python and uses NumPy for numerical computations. It simulates the motion of charged particles under the influence of electric and magnetic fields using the Boris algorithm, which is a widely used method for integrating the equations of motion for charged particles.
The results of the simulation are visualized using Matplotlib, and are validated against expected results to ensure accuracy, as well as comparison with the Runge-Kutta method. 


In nuclear fusion, the boris pusher algorithm is used to simulate the charged particles within plasma, and is incredibly effective as it conserves |u|, the magnitude of the proper velocity,  making it suitable for long-term simulations.

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

We can set this up into a system of ordinary differential equations (ODEs) that describe the motion of the particle. The Boris algorithm is then used to numerically integrate these ODEs over time, allowing us to simulate the trajectory of the charged particle in the presence of electric and magnetic fields.

$$
\begin{aligned}
   \mathbf{F} = \frac{d\mathbf{p}}{dt} = q(\mathbf{E} + \mathbf{v} \times \mathbf{B}) \\
   \frac{d\mathbf{p}}{dt} = m\frac{d\mathbf{v}}{dt}  \\
   \frac{d\mathbf{v}}{dt} = \frac{q}{m}(\mathbf{E} + \mathbf{v} \times \mathbf{B})
\end{aligned}
$$  

In this algorithm, we define the velocity of the particle as $\mathbf{u} = \gamma \mathbf{v}$, otherwise known as the "proper velocity"
Therefore the equations of motion can be rewritten in terms of the proper velocity (also known as the relativistic form of the Lorentz force) as:

$$\begin{aligned}
    \frac{d\mathbf{u}}{dt} = \frac{q}{m}(\mathbf{E} + \frac{\mathbf{u}}{\gamma} \times \mathbf{B}) \\
    \frac{d\mathbf{x}}{dt} = \frac{\mathbf{u}}{\gamma}
\end{aligned}$$


For each step over time 
$ \Delta t $
, the algorithm takes in the initial proper velocity $\mathbf{u}$ and performs the first half-kick from the electric field:

$$
\begin{aligned}
    \mathbf{u}^{-} = \mathbf{u} + \frac{q\Delta t}{2m}\mathbf{E}
\end{aligned}$$

It then performs half of the rotation due to the magnetic field:

$$
\begin{aligned}
    \mathbf{t} = \frac{q\Delta t}{2m\gamma^-}\mathbf{B} \\
    \mathbf{u}' = \mathbf{u}^{-} + \mathbf{u}^{-} \times \mathbf{t}
\end{aligned}
$$

Let us define $\mathbf{u}^{+}$ as the final proper velocity after the rotation.
To calculate the second half of the rotation, the second rotation factor $\mathbf{s}$ found by the rotation condition:

$$
\begin{aligned}
 \mathbf{u}^{+} - \mathbf{u}^{-} = (\mathbf{u}^{+} + \mathbf{u}^{-}) \times t
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
1. Calculate the half-kick from the electric field:
   $\mathbf{u}^{-} = \mathbf{u} + \frac{q\Delta t}{2m}\mathbf{E}$
2. Calculate the rotation vectors:
   $\mathbf{t} = \frac{q\Delta t}{2m\gamma^-}\mathbf{B}$
   $\mathbf{s} = \frac{2\mathbf{t}}{1+\mathbf{t}^2}$
3. Calculate the first half of the rotation:
   $\mathbf{u}' = \mathbf{u}^{-} + \mathbf{u}^{-} \times \mathbf{t}$
4. Calculate the second half of the rotation:
   $\mathbf{u}^{+} = \mathbf{u}^{-} + \mathbf{u}' \times \mathbf{s}$
5. Calculate the second half-kick from the electric field:
   $\mathbf{u}_{new} = \mathbf{u}^{+} + \frac{q \Delta t}{2m}\mathbf{E}$

Which we can update the position of the particle as:

$$
\begin{aligned}
    \mathbf{x}_{new} = \mathbf{x} + \frac{\mathbf{u}_{new}}{\gamma_{new}}\Delta t
\end{aligned}
$$
