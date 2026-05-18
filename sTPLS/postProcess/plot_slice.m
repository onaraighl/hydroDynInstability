function [X,Z,U,V,W]=plot_slice(filename)

nx=240;
nz=120;

T = readmatrix(filename, 'NumHeaderLines', 2);

X = T(:,1);
Z = T(:,2);
U = T(:,3);
V = T(:,4);
W = T(:,5);
% P = T(:,6);


X=reshape(X,nx,nz);
Z=reshape(Z,nx,nz);
U=reshape(U,nx,nz);
V=reshape(V,nx,nz);
W=reshape(W,nx,nz);

end