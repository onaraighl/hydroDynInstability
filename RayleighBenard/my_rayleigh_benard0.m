function Ra=my_rayleigh_benard0(k,Ra_guess)

Ra=fzero(@myfun,Ra_guess);

function y=myfun(x)

    tau=(x/k^4)^(1/3);
    q0=k*sqrt(tau-1);
    temp=sqrt(1+tau+tau^2)+(1+(1/2)*tau);
    q1=k*sqrt(temp/2);
    temp=sqrt(1+tau+tau^2)-(1+(1/2)*tau);
    q2=k*sqrt(temp/2);
    
    q=q1+sqrt(-1)*q2;
    
    y=imag(  (sqrt(3)+sqrt(-1))*q*tanh(q/2) )+q0*tan(q0/2);
end

end