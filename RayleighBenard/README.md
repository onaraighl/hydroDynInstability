# RayleighBenard



This directory contains some simple codes to evaluate the dispersion relation in Rayleigh-Benard convection.



# my_rayleigh_benard0



This code computes the critical Rayleigh number for the onset of instability, for a given wavenumber k (even case).  An initial guess for the critical Rayleigh number needs to be supplied.  The code corresponds to <b>Algorithm 1.1</b> in the reference text.



# my_rayleigh_benard1



This code generates the neutral curve for the onset of instability (even case).  Inputs: null.  Outputs: k\_vec and Ra\_vec.  Here, k\_vec refers to an array of wavenumbers, and Ra\_vec is the corresponding array of critical Rayeligh numbers.  The curve can be visualized by plotting:



`plot(k_vec,Ra_vec)`

This code corresponds to <b>Algorithm 1.2</b> in the reference text.

# get_eval

This code is a root-finding algorithm to obtain the allowed value of sigma (=sigma_eig) for the Rayleigh-Benard problem -- Even case.  


The code takes in the prescirbed values of Rayleigh number (Ra), Prandtl number (Pr), and wavenumber (k).  It also requires an initial guess for sigma_eig.

The code then sets up the cubic polynomial equation to solve and the resulting determinant equation.  The determinant equation is to be solved to give the allowed value of sigma_eig.  The determinant equation is solved using root-finding.


