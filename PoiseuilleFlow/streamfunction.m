function [z,phi]=streamfunction(L_aug,M_aug)
   
[~,N_1]=fix_all_parameters();

kk_1=1:(N_1-3);
eta_1=cos((pi/(N_1-2))*kk_1);
z=(1/2)*(eta_1+1);

%**************************************************************************

[V,lambda_mx]=eig(L_aug,M_aug);

lambda_eig=0*length(lambda_mx);

for i=1:length(lambda_mx)
    lambda_eig(i)=lambda_mx(i,i);
end

for i=1:length(lambda_mx)
    if(real(lambda_eig(i))>10000)
        lambda_eig(i)=-10000;
    end
end

[~,ix]=max(real(lambda_eig));
Veig=V(:,ix);


phi=0*(1:length(z));

for j=1:length(z)
   
   sum_var0=0;
   
   for k=1:(N_1+1)
       sum_var0=sum_var0+Veig(k)*T(eta_1(j),k-1);
   end
   
   phi(j)=sum_var0;
end


end

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
%**************************************************************************

