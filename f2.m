function y=f2(n,q,x)
% the two-dimensional convection-diffusion equation
%-(u_xx+u_yy)+q(u_x+u_y)+exp(u)=0.

h=1/(n+1);
r=q*h/2;
t1=2;
t2=-1-r;
t3=-1+r;
Tx=diag(t1*ones(n,1),0)+diag(t2*ones(n-1,1),-1)+diag(t3*ones(n-1,1),+1);
Ty=diag(t1*ones(n,1),0)+diag(t2*ones(n-1,1),-1)+diag(t3*ones(n-1,1),+1);


A=kron(Tx,eye(n))+kron(eye(n),Tx);

y=A*x+h*h*exp(x);
