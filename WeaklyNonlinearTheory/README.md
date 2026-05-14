# AbsoluteInstability

In <b>Chapter 7</b> of the reference text, various 1D partial differential equations (PDEs) with periodic boundary conditions are investigated.  The equations are inspired by applications in Physics, meaning they have a well-defined threshold for the onset of linear instability.  Based on the PDEs, the Fourier amplitudes (or normal modes) of the solution are computed, and are found to solve a system of coupled ODEs.  Typically, there are infinitely many degrees of freedom in the system of coupled ODEs.  However, just beyond the threshold for the onset of linear instability, an approximation can be made which greatly simplifies this coupled system.  This approxiation goes by the name of <b>weakly nonlinear theory</b>.

In this repository, the Cahn-Hilliard equation in one spatial dimension is introduced as a use case of weakly nonlinear analysis.  The equation reads:

$$
\frac{\partial C}{\partial t}=D\frac{\partial^2}{\partial x^2}(C^3-C-\gamma \partial_{xx}C),\qquad t>0,\qquad x\in (0,L),
$$

with initial data $C(x,t=0)=C_{init}(x)$ and periodic boundary coditions on the interval $(0,L)$.  Also, $D$ and $\gamma$ are positive constants.

The solution $C(x,t)$can be written in term sof Fourier modes:

$$
C(x,t)=\sum_{n=-\infty}^\infty A_n(t)e^{i k_n x},\qquad k_n=(2\pi/L)n,\qquad A_{-n}=A_n^*
$$

hence

$$
A_n=\frac{1}{L}\int_0^L C(x,t)e^{-i k_n x}dx
$$

In terms of Fourier amplitudes, the Cahn-Hilliard equation becomes:

$$
\frac{dA_n}{dt}=\nu(k_n)A_n-Dk_n^2\sum_{p=-\infty}^\infty \sum_{q=-\infty}^\infty A_p A_q A_{n-p-q},
$$

where $\nu(k)=Dk^2(1-\gamma k^2)$ gives the dispersion relation for the linearized Cahn-Hilliard equation.
