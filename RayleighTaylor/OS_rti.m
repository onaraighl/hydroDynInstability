function [lambda,L_aug,M_aug] =OS_rti(alpha)

%  28/04/2008 Solves the OS equation for the RTI

global m

global r_bottom
global r_top

global nuL
global nuG
global St
global Grav
global ReG
global im
global N_1
global N_21
global dedz_1
global dedz_21

%  The complex unit
im=sqrt(-1);

nuL=2;
nuG=2;

m=1;

r_bottom=1;
r_top=100;

St=0;
Grav=1;
ReG=1000;

N_1=100;
N_21=100;

%  Collocation grids

kk_1=1:(N_1-3);
eta_1=cos((pi/(N_1-2))*kk_1);

kk_21=1:(N_21-3);
eta_21=cos((pi/(N_21-2))*kk_21);

dedz_1 = 2/nuL;
dedz_21 = 2/nuG;

%**************************************************************************
% derive Chebyshev collocation matrices for the liquid

sizezero=zeros(N_1-3,N_1+1);
T4_1=sizezero;
T2_1=sizezero;
T0_1=sizezero;

for i=1:(N_1-3)
    for k=0:N_1         
        T4_1(i,k+1)=(dedz_1^4)*d4T(eta_1(i),k);
        T2_1(i,k+1)=(dedz_1^2)*d2T(eta_1(i),k);
        T0_1(i,k+1)=T(eta_1(i),k);          
    end
end

% derive Chebyshev collocation matrices for the gas

sizezero=zeros(N_21-3,N_21+1);
T4_21=sizezero;
T2_21=sizezero;
T0_21=sizezero;

for i=1:(N_21-3)
    for k=0:N_21              
        T4_21(i,k+1)=(dedz_21^4)*d4T(eta_21(i),k);
        T2_21(i,k+1)=(dedz_21^2)*d2T(eta_21(i),k);
        T0_21(i,k+1)=T(eta_21(i),k);                 
    end
end

%**************************************************************************
%  Define the eigenvalue problem

L_1=  (m/(r_bottom*ReG))*(T4_1 -2*alpha*alpha*T2_1 + (alpha^4)*T0_1);
M_1 = T2_1 - alpha*alpha*T0_1;

L_21= (1/(r_top*ReG))*    (T4_21-2*alpha*alpha*T2_21+ (alpha^4)*T0_21);
M_21 = T2_21 - alpha*alpha*T0_21;

%**************************************************************************

[L_aug,M_aug] = make_augmented_matrices(L_1,M_1,L_21,M_21,alpha);

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

function [L_aug,M_aug] = make_augmented_matrices(L_1,M_1,L_21,M_21,alpha)

global St
global Grav
global ReG
global im
global N_1
global N_21
global dedz_1
global dedz_21
global m

global r_bottom
global r_top


N_aug = (N_1+1) + (N_21+1)+1;

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

%**************************************************************************
% Liquid - gas interfacial conditions

for j=1:(N_1+1)
    %  IC 1
    L_aug(N_1,j)=T(1,j-1);
    
    %  IC 2
    L_aug(N_1+1,j)=dedz_1*d1T(1,j-1);  
   
    %  IC 3
    L_aug(N_1+2,j)=m*((dedz_1^2)*d2T(1,j-1)+alpha*alpha*T(1,j-1));

    %  IC 4
    L_aug(N_1+3,j)=m*((dedz_1^3)*d3T(1,j-1) - 3*alpha*alpha*dedz_1*d1T(1,j-1));
    M_aug(N_1+3,j)=r_bottom*ReG*dedz_1*d1T(1,j-1);  
end

for j=1:(N_21+1)
    %  IC 1
    L_aug(N_1,N_1+1+j)=-T(-1,j-1);
    
    %  IC 2
    L_aug(N_1+1,N_1+1+j)= -dedz_21*d1T(-1,j-1);
    
    %  IC 3
    L_aug(N_1+2,N_1+1+j)=-((dedz_21^2)*d2T(-1,j-1)+alpha*alpha*T(-1,j-1));
    
    %  IC 4
    L_aug(N_1+3,N_1+1+j)= -((dedz_21^3)*d3T(-1,j-1) - 3*alpha*alpha*dedz_21*d1T(-1,j-1));
    M_aug(N_1+3,N_1+1+j)= -ReG*r_top*dedz_21*d1T(-1,j-1);      
end

L_aug(N_1+3,N_aug)=im*alpha*ReG*((r_top-r_bottom)*Grav-alpha*alpha*St);

%**************************************************************************
% Gas block

for i=1:(N_21-3)
    for j=1:(N_21+1)
        L_aug(i+N_1+3,j+N_1+1)=L_21(i,j);
        M_aug(i+N_1+3,j+N_1+1)=M_21(i,j);
    end
end

%**************************************************************************

% Wall-gas BCs

for j=1:(N_21+1)
   L_aug(N_1+N_21+1,N_1+1+j)=T(1,j-1);
   L_aug(N_1+N_21+2,N_1+1+j)=d1T(1,j-1);
end

% Interfacial condition

for j=1:(N_1+1)
    L_aug(N_aug,j)=-im*alpha*T(1,j-1);
end

M_aug(N_aug,N_aug)=1;


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

