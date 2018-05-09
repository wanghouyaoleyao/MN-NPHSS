function DA=df2(n,q,x)


h=1/(n+1);
r=q*h/2;
t1=2;
t2=-1-r;
t3=-1+r;
Tx=diag(t1*ones(n,1),0)+diag(t2*ones(n-1,1),-1)+diag(t3*ones(n-1,1),+1);
Ty=diag(t1*zeros(n,1),0)+diag(t2*ones(n-1,1),-1)+diag(t3*ones(n-1,1),+1);
A=kron(Tx,eye(n))+kron(eye(n),Tx);

DA=A+h*h*diag(exp(x));
