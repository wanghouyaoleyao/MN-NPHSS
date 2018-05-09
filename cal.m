n=40;
q=0.2;
h=1/(n+1);
r=q*h/2;
t1=2;
t2=-1-r;
t3=-1+r;
Tx=diag(t1*ones(n,1),0)+diag(t2*ones(n-1,1),-1)+diag(t3*ones(n-1,1),+1);
Ty=diag(t1*ones(n,1),0)+diag(t2*ones(n-1,1),-1)+diag(t3*ones(n-1,1),+1);


A=kron(Tx,eye(n))+kron(eye(n),Ty);
H=(A+A')/2;
ma=max(eig(H))
mi=min(eig(H))
alpha=sqrt(ma*mi)