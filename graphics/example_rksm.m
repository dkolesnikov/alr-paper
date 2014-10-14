nh=30;
T=diag(2*ones(nh,1))+diag(-ones(nh-1,1),1)+diag(-ones(nh-1,1),-1);
I=speye(nh);
A=-(kron(T,I)+kron(I,T));
n=nh^2;
E=spdiags(rand(n,1),0,n,n);
LE=chol(E,'lower');
B=randn(n,2);
m=100;
tol=1e-9;
%tolY=1e-12;
tolY = 0.0;
opts.tol=1e-2;
s1=eigs(-A,E,1,'lm',opts);
s2=eigs(-A,E,1,'sm',opts);
ch=0;
[Z, resnorm, Zall] = rksm(A,E,LE,B,m,tol,s1,s2,ch,tolY);

fprintf('final true absolute residual norm: \n')
%disp(norm(A*Z*Z'*E+E*Z*Z'*A'+B*B')/(2 * norm(A, 'fro') * norm(Z * Z', 'fro') + norm(B * B', 'fro') ))    %this matrix should never be formed for n large
