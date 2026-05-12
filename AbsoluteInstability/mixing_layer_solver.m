function [lambda,L_aug,M_aug] =mixing_layer_solver(alpha,T4_1,T2_1,T0_1,T0_U2_1,T2_U0_1,T0_U0_1)

%  30/01/2012 Solves the OS equation for a mixing layer.

global im
global height_global
global R_param_global
global N_param_global
global Re_global
global N_global

[height,R_param,N_param,Re,N]=fix_all_parameters();

height_global=height;
R_param_global=R_param;
N_param_global=N_param;
Re_global=Re;
N_global=N;

im=sqrt(-1);

%**************************************************************************
%  Define the eigenvalue problem

VV=T4_1-2*alpha*alpha*T2_1+(alpha^4)*T0_1;

L_1=(1/Re)*VV  +   im*alpha*T0_U2_1 - im*alpha*(T2_U0_1-alpha*alpha*T0_U0_1);
M_1 = T2_1 - alpha*alpha*T0_1;

%**************************************************************************

[L_aug,M_aug] = make_augmented_matrices(L_1,M_1);

lambda=eig(L_aug,M_aug);

for i=1:length(lambda)
    if(abs(lambda(i))>100)
        lambda(i)=-10000;
    end
    
    B=isnan(lambda(i));
    if(B==1)
         lambda(i)=-10000;
    end

end



%**************************************************************************

end

%**************************************************************************
%**************************************************************************

function [L_aug,M_aug] = make_augmented_matrices(L_1,M_1)

global N_global

N=N_global;

N_aug = (N+1);

L_aug = zeros(N_aug,N_aug);
M_aug = zeros(N_aug,N_aug);

%**************************************************************************
% Liquid block

for j=1:N+1
   %  BC 1
   % dphi/dz=0 at z=0 (symmetry line - sinuous mode)
   L_aug(1,j)=d1T(-1,j-1); 
   % d3phi/dz3=0 at z=0 (va - sinuous mode)
   L_aug(2,j)=d3T(-1,j-1);
end

for i=1:(N-3)
    for j=1:(N+1)
        L_aug(i+2,j)=L_1(i,j);
        M_aug(i+2,j)=M_1(i,j);
    end
end

for j=1:N+1
    L_aug(N,j)=  T(1,j-1);
    L_aug(N+1,j)=d1T(1,j-1);
end
 


end

%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************

function t = T(z,k)

t=cos(k*acos(z));

end

%**************************************************************************

function dt = d1T(z,k)

% specifiy dt at boundary

if(abs(z)==1)
    if(z==1)
        dt=k*k;
    else
        dt=((-1)^(k+1))*k*k;
    end
else
    dt = (k/sqrt(1-z*z))*sin(k*acos(z));
end

end

%**************************************************************************

function d2t = d2T(z,k)

t=cos(k*acos(z));

if(abs(z)==1)
    if(z==1)
        d2t=(1/3)*(k^4-k^2);
    else
        d2t=(1/3)*((-1)^k)*(k^4-k^2);
    end
else
    dt = (k/sqrt(1-z*z))*sin(k*acos(z));
    d2t=(z/(1-z*z))*dt - (k*k/(1-z*z))*t;
end

end

%**************************************************************************

function d3t = d3T(z,k)

t=cos(k*acos(z));

if(abs(z)==1)
    if(z==1)
        d3t=(1/(3*5))*k*k*(k*k-1)*(k*k-2*2);
    else
        d3t=(1/(3*5))*((-1)^(k+3))*k*k*(k*k-1)*(k*k-2*2);
    end
else
    dt = (k/sqrt(1-z*z))*sin(k*acos(z));
    d2t=(z/(1-z*z))*dt - (k*k/(1-z*z))*t;
    d3t = (3*z/(1-z*z))*d2t + ((1-k*k)/(1-z*z))*dt;
end

end

%**************************************************************************

function d4t = d4T(z,k)

t=cos(k*acos(z)); 
dt = (k/sqrt(1-z*z))*sin(k*acos(z));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
d2t=(z/(1-z*z))*dt - (k*k/(1-z*z))*t;
d3t = (3*z/(1-z*z))*d2t + ((1-k*k)/(1-z*z))*dt;

d4t = (5*z/(1-z*z))*d3t + ((4-k*k)/(1-z*z))*d2t;

end

%**************************************************************************
%**************************************************************************





