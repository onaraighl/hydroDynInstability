# Spectral

This directory contains some simple codes to obtain the eigenvalues and eigenfunctions of the following equation:

$$
\frac{d^2f}{dy^2}=-\lambda f.
$$

The domain is given as $y\in (-L/2,L/2)$, and homogeneous boundary conditions are assumed on either end:

$$
f(-L/2)=f(L/2)=0.
$$


# simple

This code recasts the differential equation as ageneralized eigenvalue problem:

$$
La=\lambda Ma
$$

This equation represents a disrectization and a truncation of the differential equation with $N+1$ degrees of freedom.

The code uses numerical linear algebra in Matlab to compute the first  $N+1$ eigenvalues $\lambda$.  This code corresponds to <b>Algorithm 3.1</b> in the reference text.

# make_eigenfunction_simple

This code builds on `simple`.  The code takes the maximum eigenvalue from the set

$$
\lambda_{max}=\{\lambda_1,...,\lambda_{N_1}\}
$$

and computes the correspnoding eigenvector.



