# TransientGrowth


This directory contains a simple ODE solver to solve the following complex-valued ODE:

$$
i\frac{du}{dt}=Hu+i(\mu_0 I+G)u+a 
\begin{pmatrix}
|u_1|^2 & 0 \\
0 & |u_2|^2
\end{pmatrix}
$$

Here, $u=(u_1,u_2)^T$ is a complex-valued vector in $\mathbb{C}^2$.


