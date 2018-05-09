function [it,ot,t]=n_nphss(n,q,alpha,delta)
t=cputime;
eps = 10^(-6);
it=0;
max_iter = 1000;
x=zeros(n*n,1);       %the initial points for outer iterations(default 0)
z = zeros(n*n,1);      %the initial points for inner iterations
%alpha=q/(n+1)/2;       %the parameter in the HSS method

%      n---the number of dimension
%      q---the coefficient of the one-order-derivate term
%      alpha---the parameter in the HSS method
%      delta---the error parameter in the HSS method
%      it--the number of inner iterations
%      ot -- the number of outer iterations
%      t---the total running time 
%      x -- final iterate

errf = norm(f2(n,q,x));
for ot = 1:max_iter
    A=df2(n,q,x);  
    b=-f2(n,q,x);
    [y,itt]=nphss2(n,A,z,b,alpha,delta);
    x=y+x;
    it=it+itt;
    err = norm(f2(n,q,x))/errf;
    if (err < eps)
      break;
    end
end
t=cputime-t;

%fidhss=fopen('resulthss.txt','a');
%fprintf(fidhss,'   n     q    alpha      it     ot          time \n');
%fprintf(fidhss,'%4d  %4d %8.2f %5d %4d %9.4f\n',n,q,alpha,it,ot,t);
%fclose(fidhss);