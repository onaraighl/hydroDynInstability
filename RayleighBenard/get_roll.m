function [X,Z,U_mx,W_mx,Psi_mx]=get_roll(Ra,Pr,k,sigma_guess)

L=2*pi/k;
x=0:(L/100):L;

[~,~,z,W_vec]=get_eval_with_kernel(Ra,Pr,k,sigma_guess);
W=real(W_vec);

[X,Z]=meshgrid(x,z);
U_mx=0*X;
W_mx=0*X;
Psi_mx=0*X;

dz=z(2)-z(1);
dW=0*W;
n=length(z);

dW(1:n-1)=(W(2:n)-W(1:n-1))/dz;
dW(n)=2*dW(n-1)-dW(n-2);

for i=1:length(x)
    for j=1:length(z)
        U_mx(i,j)=-(1/k)*dW(j)*sin(k*x(i));
        W_mx(i,j)=W(j)*cos(k*x(i));
        Psi_mx(i,j)=-(1/k)*W(j)*sin(k*x(i));
    end
end

contourf(X,Z,Psi_mx.')
colormap hot
colorbar
set(gca,'fontsize',18,'fontname','arial')
set(gca,'clim',[-0.25 0.25])
xlabel('x')
ylabel('z')
drawnow

end