Let us define $\mathbf{u}^{+}$ as the final proper velocity after the rotation.
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