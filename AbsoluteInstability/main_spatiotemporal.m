function [ar,ai,omega1,omega2,omega3,omega4,omega5,omega6,omega7,omega8,omega9,omega10]=main_spatiotemporal()


% *************************************************************************
% *************************************************************************

[height,R_param,N_param,~,N]=fix_all_parameters();

kk_1=1:(N-3);
eta_1=cos((pi/(N-2))*kk_1);

y=(height*eta_1+height)/2;
% y=height*eta_1;

b=asinh(1);
sinh_val=sinh(y*b);
sinh_power=sinh_val.^(2*N_param);
shape=1./(1+sinh_power);

dsinh_power=2*N_param*b*(sinh_val.^(2*N_param-1)).*cosh(y*b);

% dshape=(-1)*dsinh_power./((1+sinh_power).^2);

ddshape1=2*(dsinh_power.^2)./((1+sinh_power).^3);
ddshape2=(2*N_param-1)*b*(sinh_val.^(2*N_param-2)).*cosh(y*b).*cosh(y*b)+...
    b*(sinh_val.^(2*N_param-1)).*sinh_val;
ddshape2=ddshape2.*(-1)*2*N_param*b./((1+sinh_power).^2);

ddshape=ddshape1+ddshape2;

ddu_vec=2*R_param*ddshape;

u_vec=1-R_param+2*R_param*shape;

dedz=2/height;
% dedz = 1/height;

%**************************************************************************
% derive Chebyshev collocation matrices for the liquid

sizezero=zeros(N-3,N+1);
T4_1=sizezero;
T2_1=sizezero;
T0_1=sizezero;
T0_U2_1=sizezero;
T2_U0_1=sizezero;
T0_U0_1=sizezero;

for i=1:(N-3)
    for k=0:N       
        
        T4_1(i,k+1)=(dedz^4)*d4T(eta_1(i),k);
        T2_1(i,k+1)=(dedz^2)*d2T(eta_1(i),k);
        T0_1(i,k+1)=T(eta_1(i),k);
                 
        T0_U2_1(i,k+1)=T(eta_1(i),k)*ddu_vec(i);
        T2_U0_1(i,k+1)=(dedz^2)*d2T(eta_1(i),k)*u_vec(i);
        T0_U0_1(i,k+1)=T(eta_1(i),k)*u_vec(i);
        
    end
end

% *************************************************************************

dalphar=0.005;
dalphai=0.02;
ar=1e-4:dalphar:1.2;
ai=-1.5:dalphai:0;

display(length(ar));

sizezero=zeros(length(ar),length(ai));
omega1=sizezero;
omega2=sizezero;
omega3=sizezero;
omega4=sizezero;
omega5=sizezero;
omega6=sizezero;
omega7=sizezero;
omega8=sizezero;
omega9=sizezero;
omega10=sizezero;

for i=1:length(ar)
    for j=1:length(ai)
        alpha_param=ar(i)+sqrt(-1)*ai(j);
        

        [lambda,~,~]=mixing_layer_solver(alpha_param,T4_1,T2_1,T0_1,T0_U2_1,T2_U0_1,T0_U0_1);
        
        [~,ix]=max(real(lambda));
        v1=sqrt(-1)*lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v2=sqrt(-1)*lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v3=sqrt(-1)*lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v4=sqrt(-1)*lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v5=sqrt(-1)*lambda(ix);
        
                lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v6=sqrt(-1)*lambda(ix);
        
                lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v7=sqrt(-1)*lambda(ix);
        
                lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v8=sqrt(-1)*lambda(ix);
        
                lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v9=sqrt(-1)*lambda(ix);
        
                lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v10=sqrt(-1)*lambda(ix);
        
        omega1(i,j)=v1;
        omega2(i,j)=v2;
        omega3(i,j)=v3;
        omega4(i,j)=v4;
        omega5(i,j)=v5;
        
        omega6(i,j)=v6;
        omega7(i,j)=v7;
        omega8(i,j)=v8;
        omega9(i,j)=v9;
        omega10(i,j)=v10;
        
        
        
    end
    display(i)
end

% *************************************************************************


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




