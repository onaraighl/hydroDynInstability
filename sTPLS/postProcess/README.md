# sTPLS/postProcess

This directory contains a suite of Fortran 90 and Matlab codes to postprocess the simulation data coming from an LES simulation using sTPLS.  

# File Structure

The main thing to keep in mind is that sTPLS outputs data every 1000 iterations, in files of the form `sp3dchannel_    N.dat`  Here, N is in integer divisible by 1000.  The white spaces are deliberate - these are to ensure the files appear in an ordered fashion in the directory hosting them.  These are very large files - it is not advisable to open them.  The first few lines are as follows:

![cartoon1.png]()

These columns correspond to (x,y,z) coordinates (first three columns), (u,v,w) velocities (next three columns).  The last two columns correspond to pressure and eddy viscosity.  

# make_slices.f90

This code reads in the files `sp3dchannel_    N.dat`.  For each N, two corresponding slices are generated: one slice in the $xz$ plane at $y=L_y/2$ and one slice in the $xy$ mindplane at $z=L_z/2$.  The slices are stored in new `.dat` files.  Two such `.dat` files are provided in this repository:

*
*

These can be visualized - e.g. the Matlab code `.m' extrasts teh xxx slice, which can then be visualized, as in the belwo cartoon.

The fortran code can be compiled in the usual way, e.g. `gfortran -O3 make_slices.f90 -o slice.x` and can then be executed by typing `/slice.x` at the command lilne.
