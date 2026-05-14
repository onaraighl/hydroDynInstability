# TransientGrowth


This directory contains a simple ODE solver to solve the following complex-valued ODE:

$$
i\frac{du}{dt}=Hu+i(\mu_0 \mathbb{I}+G)u+a 
\begin{pmatrix}
|u_1|^2 & 0 \\
0 & |u_2|^2
\end{pmatrix}
$$

Here, $u=(u_1,u_2)^T$ is a complex-valued vector in $\mathbb{C}^2$, and $a$ is a real parameter.  In this context, $H$ and $G$ are real-valued $2x2$ Hermitian matrices:

$$
H=\begin{pmatrix}
E_0 & A \\
A & E_0
\end{pmatrix}
$$

and

$$
G=\begin{pmatrix}
-g_1 & 0 \\
0 & -g_2
\end{pmatrix}
$$



