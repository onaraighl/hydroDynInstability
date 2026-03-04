function [k_vec,Ra_vec]=my_rayleigh_benard1()

k_vec=0.1:0.01:10;
Ra_vec=0*k_vec;

% For the first value of k_vec, an initial guess for the
% root is supplied by hand.  The code is quite sensitive to
% the value of this guess.

Ra=my_rayleigh_benard0(k_vec(1),500000);
Ra_vec(1)=Ra;

for i=2:length(k_vec)
    Ra_guess=Ra_vec(i-1);
    Ra=my_rayleigh_benard0(k_vec(i),Ra_guess);
    Ra_vec(i)=Ra;
end

end

    

