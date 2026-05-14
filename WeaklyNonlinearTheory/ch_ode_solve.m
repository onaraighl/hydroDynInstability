function [t_out,a1_ode,a3_ode,a5_ode]=ch_ode_solve(t_final)

% *************************************************************************
% Numerical parameters

[L,gamma,Nx,dt,n_max,my_phases,ampl_init,ampl_pert]=fix_all_parameters();

% Fundamental wavenumber:
k0=2*pi/L;

% Grid spacing:
dx=L/Nx;

% Spatial grid:
xx=(0:(Nx-1))*dx;

% % Wavenumber grid:
% kx=(2*pi/L)*(-(Nx/2):1:((Nx/2)-1));

% Final timestep:
t_ctr_final=floor(t_final/dt);

% Vector of timesteps:
t_out=0:0.05:t_final;

% *************************************************************************
% Linear Dispersion Relation:

k1=1*k0;
k3=3*k0;
k5=5*k0;

nu1=k1*k1*(1-gamma*k1*k1);
nu3=k3*k3*(1-gamma*k3*k3);
nu5=k5*k5*(1-gamma*k5*k5);

% *************************************************************************
% Initial conditions

% The initial conditions are made of a fundamental plus overtones.  The
% overtones are stored in the array pert.  Pert contains n_max overtones.
% I first of all compute the initial condition in real space, as in
% ch_oned_solve.m

pert=0*xx;

for n=2:n_max
    pert=pert+cos(n*k0*xx+2*pi*my_phases(n));
end

pert=pert/(n_max-2+1);

c_init=ampl_init*( cos(k0*xx) + ampl_pert*pert );
c_hat=fftshift(fft(c_init));

% Now I obtain the proper initial conditions in Fourier space:

% a1_wnl=(dx/L)*sum(c_init.*exp(-sqrt(-1)*k0*xx)), which explicitly
% integrates to a1_wnl=ampl_init*0.5;, since
% (dx/L)*sum(ccos(k0*xx).*exp(-sqrt(-1)*k0*xx))=0.5.

a1_wnl=ampl_init*0.5;

% *************************************************************************

[T,a1_wnlT]=ode45(@myfun, [0 t_final], a1_wnl );

a1_ode=interp1(T,a1_wnlT,t_out,'spline');
a3_ode=k3*k3*(a1_ode.^3)/nu3;
a5_ode=3*k5*k5*(a1_ode.*a1_ode.*a3_ode)/nu5;

% *************************************************************************

    function dAdt=myfun(~,a)
        
        a1_wnl=a;
        
        % Slaving step:
        a3_wnl=k3*k3*(a1_wnl^3)/nu3;
        a5_wnl=3*k5*k5*(a1_wnl*a1_wnl*a3_wnl)/nu5;

        % Computation of interaction terms for A1-equation:

        nl1=3*a1_wnl*conj(a1_wnl)+6*a3_wnl*conj(a3_wnl)+6*a5_wnl*conj(a5_wnl);
        nl1=-k1*k1*a1_wnl*nl1;

        nl2=conj(a1_wnl)*conj(a1_wnl)+a3_wnl*conj(a5_wnl);
        nl2=-3*k1*k1*a3_wnl*nl2;

        nl3=conj(a1_wnl)*conj(a3_wnl);
        nl3=-6*k1*k1*a5_wnl*nl3;
        
        % Expression for dA_1/dt:
        dAdt=nu1*a1_wnl+(nl1+nl2+nl3);

    end

end
