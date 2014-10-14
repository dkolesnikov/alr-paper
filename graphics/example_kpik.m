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
tolY=1e-12;
[Z,r]=kpik(A,E,LE,B,m,tol,tolY);

fprintf('final true absolute residual norm: \n')
disp(norm(A*Z*Z'*E+E*Z*Z'*A'+B*B'))    %this matrix should never be formed for n large 
disp(size(Z))
