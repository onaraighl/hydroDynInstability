function [lambda,L_aug,M_aug] =OS_solver(alpha,u_vec1,ddu_vec1)

%  02/07/2012 Solves the OS equation for a Couette flow in a channel.

global im
global Re_global
global N_1_global

[Re,N_1]=fix_all_parameters();

Re_global=Re;
N_1_global=N_1;

im=sqrt(-1);

kk_1=1:(N_1-3);
eta_1=cos((pi/(N_1-2))*kk_1);

% dedz=2 for channel flow in a channel z\in[0,1].
dedz_1 = 2;

%**************************************************************************
% derive Chebyshev collocation matrices for the liquid

sizezero=zeros(N_1-3,N_1+1);
T4_1=sizezero;
T2_1=sizezero;
T0_1=sizezero;
VV  =sizezero;
T0_U2_1=sizezero;
T2_U0_1=sizezero;
T0_U0_1=sizezero;

for i=1:(N_1-3)
    for k=0:N_1       
        
        T4_1(i,k+1)=(dedz_1^4)*d4T(eta_1(i),k);
        T2_1(i,k+1)=(dedz_1^2)*d2T(eta_1(i),k);
        T0_1(i,k+1)=T(eta_1(i),k);
         
        VV(i,k+1)=(dedz_1^4)*d4T(eta_1(i),k)-2*alpha*alpha*(dedz_1^2)*d2T(eta_1(i),k)+(alpha^4)*T(eta_1(i),k);
        
        T0_U2_1(i,k+1)=T(eta_1(i),k)*ddu_vec1(i);
        T2_U0_1(i,k+1)=(dedz_1^2)*d2T(eta_1(i),k)*u_vec1(i);
        T0_U0_1(i,k+1)=T(eta_1(i),k)*u_vec1(i);
        
    end
end


%**************************************************************************
%  Define the eigenvalue problem

L_1=(1/Re)*VV  +   im*alpha*T0_U2_1 - im*alpha*(T2_U0_1-alpha*alpha*T0_U0_1);
M_1 = T2_1 - alpha*alpha*T0_1;

%**************************************************************************

[L_aug,M_aug] = make_augmented_matrices(L_1,M_1);

lambda=eig(L_aug,M_aug);

for i=1:length(lambda)
    if(abs(lambda(i))>10000)
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
%**************************************************************************
%**************************************************************************

function [L_aug,M_aug] = make_augmented_matrices(L_1,M_1)

global N_1_global

N_1=N_1_global;

N_aug = (N_1+1);

L_aug = zeros(N_aug,N_aug);
M_aug = zeros(N_aug,N_aug);

%**************************************************************************
% Liquid block

for j=1:N_1+1
   %  BC 1
   L_aug(1,j)=T(-1,j-1); 
   %  BC 2
   L_aug(2,j)=d1T(-1,j-1);
end

for i=1:(N_1-3)
    for j=1:(N_1+1)
        L_aug(i+2,j)=L_1(i,j);
        M_aug(i+2,j)=M_1(i,j);
    end
end

for j=1:N_1+1
   %  BC 4
   L_aug(N_1,j)=T(1,j-1);
   %  BC 4
   L_aug(N_1+1,j)=d1T(1,j-1);
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

% % function d3t = d3T(z,k)
% % 
% % t=cos(k*acos(z));
% % 
% % if(abs(z)==1)
% %     if(z==1)
% %         d3t=(1/(3*5))*k*k*(k*k-1)*(k*k-2*2);
% %     else
% %         d3t=(1/(3*5))*((-1)^(k+3))*k*k*(k*k-1)*(k*k-2*2);
% %     end
% % else
% %     dt = (k/sqrt(1-z*z))*sin(k*acos(z));
% %     d2t=(z/(1-z*z))*dt - (k*k/(1-z*z))*t;
% %     d3t = (3*z/(1-z*z))*d2t + ((1-k*k)/(1-z*z))*dt;
% % end
% % 
% % end

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


