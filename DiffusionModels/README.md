# DiffusionModels

This directory contains some simple codes to solve the <i>Model Diffusion Problem</i>, as discussed in <b>Chapter 10</b> of the reference text.

# Model Diffusion Problem

The Model Diffusion Problem refers to the following Partial Differential Equation (PDE):

$$
\frac{\partial C}{\partial t}=\nabla^2 C+s(x,z),\qquad (x,z)\in\Omega,
$$

where

$$
\Omega=(0,L_x)\times (0,1),
$$

and $\nabla^2=\partial_x^2+\partial_z^2$ is the Laplacian.

The partial differential equation is subject to the following boundary conditions:

$$
\frac{\partial C}{\partial z}=0,\qquad z=0,\qquad z=1,
$$

together with  periodic boundary conditions in the $x$ direction:

$$
C(x=0,z,t)=C(x=L_x,z,t).
$$

Finally, an initial condition is prescribed:

$$
C(x,z,t=0)=C_{\mathrm{init}}(x,z),\qquad (x,z)\in \overline{\Omega},
$$

where $C_{\mathrm{init}}(x,z)$ is a continuous function.

# Discretezation

We introduce centred differencing in space as a way of approximating the Laplace operator numerically.  We  introduce Crank-Nicolson temporal discretization as a way of discretizing the temporal derivative $\partial/\partial t$.  Crank-Nicolson is a so-called implicit method, which means that a certain equation must be inverted in order to evolve the numerical solution forward in time, stepping from one time step to the next.  

We therefore disccretize the model PDE and compute the approximate numerical solution on a discrete grid:

$$
\begin{aligned}
x_i&=(i-1)\Delta x,\qquad i=1,\cdots n_x,\\
z_j&=(j-1)\Delta z,\qquad j=1,\cdots n_z,
\end{aligned}
$$

such that

$$
(n_x-1) \Delta x=L_x, \qquad \Delta x=L_x/(n_x-1),
$$

and similarly, $\Delta z=L_z/(n_z-1)$.

The PDE is also discretized in time, such that the solution is only available at discrete points in time $t_n=n\Delta t$, with $n=0,1,\cdots$.  The solution at $t_n$ and $\boldsymbol{x}=(i\Delta x,\Delta j)$ is written as $C_{ij}^n$.  The diffusion operator in the PDE is approximated by centred differences:

$$
\begin{aligned}
\left(\nabla^2C\right)_{ij}\approx \frac{C_{i+1,j}+C_{i-1,j}-2C_{ij}}{\Delta x^2}+\frac{C_{i,j+1}+C_{i,j-1}-2C_{ij}}{\Delta z^2}:=\mathcal{D}(C_{ij})\\
i=2,3,\cdots,n_x-1,\qquad j=2,3,\cdots,n_z-1.
\end{aligned}
$$

The discretization in time is done using a Crank-Nicolson scheme:

$$
\frac{C_{ij}^{n+1}-C_{ij}^n}{\Delta t}=\tfrac{1}{2}\left[\mathcal{D}(C^n_{ij})+\mathcal{D}(C^{n+1}_{ij})\right]+s_{ij},
$$

valid for $i=2,3,\cdots,n_x-1$ and $j=2,3,\cdots,n_z-1$.  

We re-arrange the discretize diffusion equation as follows:

$$
\left[\mathbb{I}-\tfrac{1}{2}\Delta t\mathcal{D}\right]\left(C_{ij}^{n+1}\right)=\left[\mathbb{I}+\tfrac{1}{2}\Delta t\mathcal{D}\right]\left(C_{ij}^n\right)+\Delta ts_{ij}.
$$

On the left-hand side, the quantity $\left[1-\tfrac{1}{2}\Delta t\mathcal{D}\right]$ is in fact a matrix operator, and the solution is available only in <i>implicit</i> form: an inversion needs to be performed to extract $C_{ij}^{n+1}$ from this implicit equation:

$$
C_{ij}^{n+1}=\left[\mathbb{I}-\tfrac{1}{2}\Delta t\mathcal{D}\right]^{-1}\{\left[\mathbb{I}+\tfrac{1}{2}\Delta t\mathcal{D}\right]\left(C_{ij}^n\right)+\Delta ts_{ij}\}.
$$

The implicit equation is written out in more detail now:

$$
\left(1+a_x+a_z\right)C_{ij}^{n+1}-\tfrac{1}{2}a_x\left(C_{i+1,j}^{n+1}+C_{i-1,j}^{n+1}\right)-
\tfrac{1}{2}a_z\left(C_{i,j+1}^{n+1}+C_{i,j-1}^{n+1}\right)=
\left[\mathbb{I}+\tfrac{1}{2}\Delta t\mathcal{D}\right]\left(C_{ij}^n\right)+s_{ij}:=\mathrm{RHS}_{ij}^n,
$$

where $a_x=\Delta t/\Delta x^2$ and $a_z=\Delta t/\Delta z^2$.

We tidy up the above expression as follows:

$$
\left(1+a_x+a_y\right)C_{ij}^{n+1}-\tfrac{1}{2}a_x\left(C_{i+1,j}^{n+1}+C_{i-1,j}^{n+1}\right)-
\tfrac{1}{2}a_z\left(C_{i,j+1}^{n+1}+C_{i,j-1}^{n+1}\right)=\mathrm{RHS}_{ij}^n.
$$

This last expression is the one that gets implmeneted in Matlab.  The matrix inversion is achieved using either the Jacobi method or the SOR method.

# Matlab codes

The above numerical algorithm (with Jacobi for the matrix inversion) is implemented in the following Matlab codes:

* `test_diffusion_jacobi.m`
*  `fix_all_parameters.m`

The first file here is the main one.  The inputs are null, as key parameters are set in `fix_all_parameters.m`.  The outputs are:

* xx - array of discrete x-coordinates
* yy - array of discrete y-coordinates
* Concentration field C at the final time.
* time_vec - array of discrete times
* norm_decay - L2 norm of the differene between $C(x,y,t)$ and $C_0(x,y)$, the latter being the solutoin of the stationary problem.

Arrays `time_vec` and `norm_decay` are to be taken together, thus `norm_decay` can be visualized by plotting e.g. `plot(time_vec,norm_decay)`

# Stationary Problem

At late times, the solution of the model diffusion problem relaxes to the solution of a corresponding model Poisson problem - this is the already-introduced function $C_0(x,y)$.  The code `test_poisson_jacobi.m`  solves this Poisson problem.
