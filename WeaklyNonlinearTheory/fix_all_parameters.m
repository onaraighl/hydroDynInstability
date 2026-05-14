function [L,gamma,Nx,dt,n_max,my_phases,ampl_init,ampl_pert]=fix_all_parameters()

L=1;
gamma=(L*L/(8*pi*pi));

% Number of gridpoints:
Nx=256;

% Timestep:
dt=1e-5;

% Setup of initial condition:
n_max=10;

% Ten random phases:
my_phases=[ 0.1576    0.9706    0.9572    0.4854    0.8003    0.1419    0.4218    0.9157    0.7922    0.9595];
ampl_init=1e-4;
ampl_pert=1e-1;


end