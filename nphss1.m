function [x,it]= nphss1(n,A, x, b, alpha,delta)     
% nphss1.m solves the linear system Ax=b using the NPHSS Method for inner iteration!


% input   n        the number of dimension
%         A        real matrix,that is F'(x^k)
%         x        real initial guess vector,generally, 0 is the default value.
%         b        real right hand side vector,that is -F(x^k)
%         alpha    the parameter in NPHSS method
%         delta    real, error tolerance for inner iteration,0.01,0.1,0.5!    
%
% output  x        real solution vector
%         error    real error norm,when norm(b-a*x)/norm(b-a*x^0) is no greater 
%                  than \delta,the inner iteration is terminated!
%         it       integer number of inner iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it
  flag = 0;             % initialization
  max_it=1000;          %integer,maximum number of iterations 
  P=diag(diag(A));
%  bnrm2 = norm( b );
%  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end
  r = b - A*x;
  nr=norm(r);
%  error = norm( r ) / bnrm2;
%  if ( error < tol ) return, end

  H=(A+A')/2;
  S=(A-A')/2;

 [LS,US]=lu(alpha*P+H);

  for iter = 1:max_it                         % begin iteration
     x_1 = x;
     x=US\(LS\((alpha*P-S)*x+b));
     x=US\(LS\((alpha*P-S)*x+b));
     error = norm( b - A*x ) / nr;     % compute error
     if ( error <= delta ) it=iter;break, end          % check convergence
  end
  it=iter;                                
  if ( error > delta ) flag = 1; end;           % no convergence

% END hss.m