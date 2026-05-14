function [t,u1,u2,energy,t0,energy0]=my_ode_modelB(a,b,ampl_init,tfin)

% *************************************************************************

% Performs DNS of model B (nonlocal).
% 
% Input paramters are initial conditions i.t.o. eigenvalues of linear
% problem, as well as initial amplitude and final time.

% *************************************************************************


my_i=sqrt(-1);

small=0.01;

A=.1;
E0=1;

g2=0.01;
g1=g2+sqrt(4*A*A+small*small);

mu_min=min(g1,g2);
mu_max=0.5*(g1+g2)-sqrt( (g1-g2)*(g1-g2)-4*A*A);

mu0=0.95*mu_max;

H=[E0,A;
  A,E0];
G=[mu0-g1,0;0,mu0-g2];
L=H+my_i*G;

a_nonlin=0;

% *************************************************************************

t0=0:0.1:20;
energy0=0*t0;

for i=1:length(t0);
    A=exp(my_i*(L')*t0(i))*exp(-my_i*L*t0(i));
    [V,~]=eig(A);
    v1=V(:,1);
    v2=V(:,2);
    energy0(i)=sqrt(max(norm(A*v1),norm(A*v2)));
end

% *************************************************************************

init=[a,b].';
init=ampl_init*init/norm(init);


%options = odeset('RelTol',1e-16,'AbsTol',1e-20);

[T,Y]=ode87(@my_ode_loc,[0,tfin],init);

u1=Y(:,1);
u2=Y(:,2);
t=T;

energy=sqrt(u1.*conj(u1)+u2.*conj(u2));
% energy=energy/energy(1);
% u1=u1/sqrt(energy(1));
% u2=u2/sqrt(energy(1));

function dydt=my_ode_loc(~,y)
    nonlin_mx=a_nonlin*diag([y(1)*conj(y(1)),y(2)*conj(y(2))]);
    dydt=-my_i*L*y-my_i*nonlin_mx*y;
end

end