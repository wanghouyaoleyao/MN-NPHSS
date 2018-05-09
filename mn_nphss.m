
function [it1,it2,ot,t]=mn_nphss(n,q,alpha,delta)
t=cputime;
eps = 10^(-6);
it1=0;
it2=0;
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
    b1=-f2(n,q,x);
    [d,itt]=nphss2(n,A,z,b1,alpha,delta);
    y=x+d;
    it1=it1+itt;
   
    
    b2=-f2(n,q,y);
    [h,itt]=nphss2(n,A,z,b2,alpha,delta);
    x=y+h;
    it2=it2+itt;
    
    err = norm(f2(n,q,x))/errf;
    if (err < eps)
      break;
    end
end
t=cputime-t;
end
