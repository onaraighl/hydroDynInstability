function [sigma_eig1,sigma_eig2,test]=get_eval(Ra,Pr,k,sigma_guess)

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

sigma_eig1=fzero(@myfun1,sigma_guess);
sigma_eig2=fzero(@myfun2,sigma_guess);
test=abs(sigma_eig2-sigma_eig1);


%**************************************************************************
% Definition of the function myfun.

    function y1=myfun1(sigma)

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
  
        % Define the function whose zeros are to be extracted.
        % The real and imaginary parts of the determinant should vanish
        % simultaneously.  The correct condition is
        %
        % y=[real(det(M))]^2+[imag(det(M))]^2, y=0.       (1)
        %
        % Finding y=0 from Equation (1) is optimization, not rootfinding.
        % By trial and error, we have found that the optization in matlab
        % is difficult, even when using fminbnd.  We have therefore chosen
        % instead to break (1) into two problems:
        %
        % y1=real(det(M))+imag(det(M));                    (2)
        % y2=real(det(M))-imag(det(M));                    (3)
        %
        % If y1(sigma1)=0 and y1(sigma2=0) and if sigma1=sigma2, then the
        % simultaneous solution of (2) and (3) is equivalent to the
        % solution of (1).  Hence,
        
        y1=real(det(M))+imag(det(M));
        
        % We define y2 in the code below.
    end

    function y2=myfun2(sigma)

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
            
            y2=real(det(M))-imag(det(M));
    end



end