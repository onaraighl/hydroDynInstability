function [sigma_eig,M_ker,z_vec,W_vec]=get_eval_with_kernel(Ra,Pr,k,sigma_guess)

% Rootfinding solver to obtain the allowed value of sigma (=sigma_eig) for
% the Rayleigh-Benard problem -- Even case.
%
% The code takes in the prescirbed values of Rayleigh number (Ra), 
% Prandtl number (Pr), and wavenumber (k).  It also requires an initial
% guess for sigma_eig.
%
% The code then sets up the cubic polynomial equation to solve and the
% resulting determinant equtaion.  The determinant equation is to be solved
% to give the allowed value of sigma_eig.
%
% The determinant equation is solved using rootfinding.

ksq=k*k;

% Call to Matlab's fzero solver using the initial guess sigma_guess and the
% function handel myfun.

sigma_eig=fzero(@myfun,sigma_guess);

%**************************************************************************
% Setup of the matrix M once the eigenvalue sigma_eig is known - the idea
% now is that the resulting matrix has a non-trivial kernel/nullspace.

p=[1,-(sigma_eig+sigma_eig*Pr), sigma_eig*sigma_eig*Pr,Ra*k*k];
y_roots=roots(p);
x_roots=y_roots+k*k;    
q1=sqrt(x_roots(1));
q2=sqrt(x_roots(2));
q3=sqrt(x_roots(3));

M_ker(1,1)=cosh(q1/2);
M_ker(1,2)=cosh(q2/2);
M_ker(1,3)=cosh(q3/2);

M_ker(2,1)=q1*sinh(q1/2);
M_ker(2,2)=q2*sinh(q2/2);
M_ker(2,3)=q3*sinh(q3/2);

M_ker(3,1)=(q1*q1-ksq)*( (q1*q1-ksq)-sigma_eig)*cosh(q1/2);
M_ker(3,2)=(q2*q2-ksq)*( (q2*q2-ksq)-sigma_eig)*cosh(q2/2);
M_ker(3,3)=(q3*q3-ksq)*( (q3*q3-ksq)-sigma_eig)*cosh(q3/2);

% r=null(M_ker);
[V,D]=eig(M_ker);
d=diag(D);
[~,ix]=min(abs(d));
r=V(:,ix);

z_vec=-0.5:0.01:0.5;
W_vec=r(1)*cosh(q1*z_vec)+r(2)*cosh(q2*z_vec)+r(3)*cosh(q3*z_vec);

%**************************************************************************
% Definition of the function myfun.

    function y=myfun(sigma)

        % Set up the cubic polynomial.
        p=[1,-(sigma+sigma*Pr), sigma*sigma*Pr,Ra*k*k];
        % Extract polynomial roots.
        y_roots=roots(p);
        % Back out the allowed values of q:
        x_roots=y_roots+k*k;    
        q1=sqrt(x_roots(1));
        q2=sqrt(x_roots(2));
        q3=sqrt(x_roots(3));

        % Set up the matrix M, as in the book.
        M(1,1)=cosh(q1/2);
        M(1,2)=cosh(q2/2);
        M(1,3)=cosh(q3/2);
        
        M(2,1)=q1*sinh(q1/2);
        M(2,2)=q2*sinh(q2/2);
        M(2,3)=q3*sinh(q3/2);
        
        M(3,1)=(q1*q1-ksq)*( (q1*q1-ksq)-sigma)*cosh(q1/2);
        M(3,2)=(q2*q2-ksq)*( (q2*q2-ksq)-sigma)*cosh(q2/2);
        M(3,3)=(q3*q3-ksq)*( (q3*q3-ksq)-sigma)*cosh(q3/2);
  
        y=real(det(M))+imag(det(M));

    end



end