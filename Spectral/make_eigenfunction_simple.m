function [y,psi]=make_eigenfunction_simple(lambda_eig,V)

global N

%**************************************************************************

N =100;
L=2*pi;

%  Collocation grids

int_range=1:(N-1);
xx=cos((pi/(N+0))*int_range);

%**************************************************************************

[mx,ix]=min(real(lambda_eig));
lambda_eig(ix)=1000;
[mx,ix]=min(real(lambda_eig));
display(mx)
Veig=V(:,ix);

% derive Chebyshev collocation matrices 

psi=0*xx;

for i=1:(N-1)
    x=xx(i);
    sum_val=0;
    
    for k=1:(N+1)
        sum_val=sum_val+Veig(k)*T(x,k-1);
    end
    
    psi(i)=sum_val;
end

y=(L/2)*xx;

miny=min(y);
maxy=max(y);

dy=(maxy-miny)/1000;

yi=miny:dy:maxy;
psii=interp1(y,psi,yi,'spline');

y=yi;
psi=psii;



%**************************************************************************

end

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

