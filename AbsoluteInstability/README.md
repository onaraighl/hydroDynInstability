# AbsoluteInstability



This directory contains a Chebyshev collocation method to compute the eigenvalues of the Orr-Sommerfeld equation in case of a mixing-layer flow, given by:

$$
U_0(z)=1-\Lambda+\frac{2\Lambda}{1+\sinh^{2N}[z\sinh^{-1}(1)]},
$$

Nominally, $z\in (-\infty,\infty)$, which is truncated in the code, with $z\in (-H,H)$, where $H$ should be large.  This example is covered in Section 5.7 of the reference text.

Again, the code structure is a little bit complicated, and is summarized here in what follows.

Calls to the core solver `mixing_layer_solver.m` are handled by various wrap-around functions, supported by:

* `fix_all_parameters.m`
* `main_spatiotemporal.m`

# fix_all_parameters

In this code, the Reynolds number $Re$ and the truncation number $N_1$ are fixed.  Correspondingly, there are $N_1+1$ Chebyshev polynomials used in the truncation.  Other parameters are also set in a self-explanatory way:

```matlab
height=8;
R_param=-1.1;
N_param=5;
Re=100;

%  Order of Chebyshev approximation

N_1 =95;
```

Here, `R_param` is equivlanet to $\Lambda$, and `N_param` is equivalent to $N$ in the definition of the base flow $U_0(z)$.

# main_spatiotemporal

In this code, a grid of complex-valued wavenumbers is created:

```matlab
dalphar=0.005;
dalphai=0.02;
ar=1e-4:dalphar:1.2;
ai=-1.5:dalphai:0;
```

A `for` loop over all such complex-valued wavenumbers is created, and each each wavenumber, a call is made to `mixing_layer_solver.m`  This code recasts the Orr-Sommerfeld problem for $U_0(z)$ as a generalized eigenvalue problem $La=\lambda Ma$, and solves for the first $N_1+1$ eigenvalues:

```matlab
[lambda,~,~]=mixing_layer_solver(alpha_param,T4_1,T2_1,T0_1,T0_U2_1,T2_U0_1,T0_U0_1);
```

The eigenvalues are sorted by the largest real part as in previous examples, and arrays are built up the complex-valued frequencies, which depend on wavenumber:

```matlab
for i=1:length(ar)
    for j=1:length(ai)
        alpha_param=ar(i)+sqrt(-1)*ai(j);        

        [lambda,~,~]=mixing_layer_solver(alpha_param,T4_1,T2_1,T0_1,T0_U2_1,T2_U0_1,T0_U0_1);
        
        [~,ix]=max(real(lambda));
        v1=sqrt(-1)*lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v2=sqrt(-1)*lambda(ix);
        lambda(ix)=-1000;
	    …
 	    omega1(i,j)=v1;
        omega2(i,j)=v2;
	    …
	end
end
```

Results can be visualized using contour plots:

`contour(ar,ai,imag(omega1))`,

etc.






