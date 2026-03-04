function [sigma_guess,res_vec]=get_eval0(Ra,Pr,k)

ksq=k*k;
% sigma_eig=fzero(@myfun,sigma_guess);

sigma_guess=-30:0.1:30;
res_vec=0*sigma_guess;

for i=1:length(sigma_guess)
    y_val=myfun(sigma_guess(i));
    res_vec(i)=y_val;
end


    function y=myfun(sigma)
        
        p=[1,-(sigma+sigma*Pr), sigma*sigma*Pr,Ra*k*k];
        y_roots=roots(p);
        x_roots=y_roots+k*k;
        
        q1=sqrt(x_roots(1));
        q2=sqrt(x_roots(2));
        q3=sqrt(x_roots(3));

        M(1,1)=cosh(q1/2);
        M(1,2)=cosh(q2/2);
        M(1,3)=cosh(q3/2);
        
        M(2,1)=q1*sinh(q1/2);
        M(2,2)=q2*sinh(q2/2);
        M(2,3)=q3*sinh(q3/2);
        
        M(3,1)=(q1*q1-ksq)*( (q1*q1-ksq)-sigma)*cosh(q1/2);
        M(3,2)=(q2*q2-ksq)*( (q2*q2-ksq)-sigma)*cosh(q2/2);
        M(3,3)=(q3*q3-ksq)*( (q3*q3-ksq)-sigma)*cosh(q3/2);
  
        y=real(det(M))^2+imag(det(M))^2;
    end



end