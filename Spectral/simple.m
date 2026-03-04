function [ll,lambda_eig,V]=simple()

global N

%**************************************************************************

N =100;
L=2*pi;

%  Collocation grids

int_range=1:(N-1);
xx=cos((pi/(N+0))*int_range);

%**************************************************************************

% derive Chebyshev collocation matrices 

size_zero=zeros(N-1,N+1);

T2=size_zero;
M_mx=size_zero;

for i=1:(N-1)
    for k=0:N
        x=xx(i);
        dy_x=2/L;
        T2(i,k+1)=(dy_x^2)*d2T(x,k);        
        M_mx(i,k+1)=-T(x,k);
    end
end

%**************************************************************************
%  Define the eigenvalue problem

L_1=T2;
M_1 = M_mx;

%**************************************************************************

[L_aug,M_aug] = make_augmented_matrices(L_1,M_1);

[V,d]=eig(L_aug,M_aug);
[nd,~]=size(d);

lambda_eig=0*(1:nd);
for i=1:nd
    lambda_eig(i)=real(d(i,i));
end

for i=1:nd
    if(abs(lambda_eig(i))>1e6);
        lambda_eig(i)=1000;
    end
end

ll=sort(lambda_eig);

%**************************************************************************

end


%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************

function [L_aug,M_aug] = make_augmented_matrices(L_1,M_1)

global N

N_aug = (N+1);

L_aug = zeros(N_aug,N_aug);
M_aug = zeros(N_aug,N_aug);

% Liquid block

for j=1:N+1
   %  BC 1
   L_aug(1,j)=T(-1,j-1); 
   M_aug(1,j)=0; 
end

for i=1:(N-1)
    for j=1:(N+1)
        L_aug(i+1,j)=L_1(i,j);
        M_aug(i+1,j)=M_1(i,j);
    end
end

for j=1:N+1
   %  BC 2
   L_aug(N+1,j)=T(1,j-1); 
   M_aug(N+1,j)=0; 
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
%**************************************************************************

