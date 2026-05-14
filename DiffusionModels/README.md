# DiffusionModels

This directory contains some simple codes to solve the <i>Model Diffusion Problem</i>, as discussed in <b>Chapter 10</b> of the reference text.

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
