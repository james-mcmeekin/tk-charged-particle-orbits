Let us split the velocity of the particle into two components:

$$
\begin{aligned}
    \mathbf{v} = \mathbf{v}_{d} + \mathbf{w}
\end{aligned}
$$

Where $\mathbf{v}_{d}$ is the drift velocity, and that the electric force and magnetic force of the drift cancel; that is:

$$
\begin{aligned}
    \mathbf{E} + \mathbf{v}_{d} \times \mathbf{B} = 0
\end{aligned}
$$

All we need to do is cross this condition with $\mathbf{B}$:

$$
\begin{aligned}
    \mathbf{E} \times \mathbf{B} + (\mathbf{v}_{d} \times \mathbf{B})\times \mathbf{B} = 0 \\
    \mathbf{E} \times \mathbf{B} + \mathbf{B}(\mathbf{v}_{d} \cdot \mathbf{B}) - \mathbf{v}_{d}\mathbf{B}^2 = 0 \\
    \mathbf{E} \times \mathbf{B} = \mathbf{v}_{d}\mathbf{B}^2 - \mathbf{B}(\mathbf{v}_{d} \cdot \mathbf{B})
\end{aligned}
$$

Next, we can split up $\mathbf{v}_{d}$ into perpendicular and parallel components $\mathbf{v}_{d}=\mathbf{v}_{d,\perp}+\mathbf{v}_{d,||}$. We note that $\mathbf{B}\cdot \mathbf{v}_{d,\perp}=0$. Substituting:

$$
\begin{aligned}
    \mathbf{E} \times \mathbf{B} = \mathbf{v}_{d}\mathbf{B}^2 - \mathbf{B}(\mathbf{v}_{d,\perp}+\mathbf{v}_{d,||} \cdot \mathbf{B}) \\
    \mathbf{E} \times \mathbf{B} = \mathbf{B}^2(\mathbf{v}_{d}-\mathbf{v}_{d,||}) \\
    \mathbf{E} \times \mathbf{B} = \mathbf{B}^2\mathbf{v}_{d,\perp} \\
    \mathbf{v}_{d,\perp} = \frac{\mathbf{E} \times \mathbf{B}}{\mathbf{B}^2}
\end{aligned}
$$