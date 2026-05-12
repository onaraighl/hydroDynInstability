function [alpha,lambda1,lambda2]=main_temporal()
% function [lambda,L,M]=main_temporal(alpha_param)

[u_vec1,ddu_vec1,~]=get_u();


% [lambda,L,M]=OS_solver(alpha_param,u_vec1,ddu_vec1);


% **************************************************************************

alpha=0:0.0001:0.005;
display(length(alpha));

sizezero=0*(1:length(alpha));
lambda1=sizezero;
lambda2=sizezero;
lambda3=sizezero;
lambda4=sizezero;
lambda5=sizezero;
lambda6=sizezero;
lambda7=sizezero;
lambda8=sizezero;
lambda9=sizezero;
lambda10=sizezero;

for i=1:length(alpha)
        alpha_param=alpha(i);
        [lambda,~,~]=OS_solver(alpha_param,u_vec1,ddu_vec1);
        [~,ix]=max(real(lambda));
        v1=lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v2=lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v3=lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v4=lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v5=lambda(ix);
        lambda(ix)=-1000;

        [~,ix]=max(real(lambda));
        v6=lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v7=lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v8=lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v9=lambda(ix);
        lambda(ix)=-1000;
        [~,ix]=max(real(lambda));
        v10=lambda(ix);
        
        lambda1(i)=v1;
        lambda2(i)=v2;
        lambda3(i)=v3;
        lambda4(i)=v4;
        lambda5(i)=v5;
        lambda6(i)=v6;
        lambda7(i)=v7;
        lambda8(i)=v8;
        lambda9(i)=v9;
        lambda10(i)=v10;
        
        display(i)
end

end




