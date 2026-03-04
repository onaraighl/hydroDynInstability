function [alpha_vec,lambda_vec]=call_OS_rti()

alpha_vec=1e-6:0.5:10;
lambda_vec=0*alpha_vec;

for i=1:length(alpha_vec)
    [lambda] =OS_rti(alpha_vec(i));
    
    [~,ix]=max(real(lambda));
    lambda_vec(i)=lambda(ix);
    display(i)
end
    


end