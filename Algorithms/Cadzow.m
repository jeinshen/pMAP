function [x, It, Error] = Cadzow( x, pars)
% This is an implementation of Cadzow(alpha) method in the paper
% Iteratiove algorithms for weighted and unweighted finite-rank time sereis
% approximation, Statistics and Its Interface, 10 (2017), 5--18
% by N. Zvonarev and N.  Golyandina
% 
% Input: 
% alpha: 0 < alpha < 1 (e.g., 0.2)
%  x -- noisy signal
%
%  pars
%     compulsary pars
%       pars.L -- window length
%       pars.r -- rank order
%     
%     optional pars
%       pars.maxit -- maximum number of iterations
%       pars.tol   -- tolerance to terminate the algorithm (1.0e-4, used by [ZG17])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start: Input the parameters
 L = pars.L;
 r = pars.r;
 
 if isfield(pars, 'maxit')
     maxit = pars.maxit;
 else
     maxit = 100; % default used by ZG17 paper
 end
 
 if isfield(pars, 'tol')
     tol = pars.tol;
 else
     tol = 1.0e-3;
 end
 
 N = length(x);
 K = N-L+1;

 % Initial point X0
  X0 = hankel(x(1:L), x(L:end));
  A = X0;
  x0 = x;
  
 % Initial error = 1
  Error = 1;
  It = 0;
  err = 1;
 % get the weight matrix W for HankelProjection.m 

  while (It <= maxit) && (Error > tol)
      M = LowRankProj(X0, r);
      X1 = HankelMatrix(M);
      
      x1 = [X1(:, 1); X1(L, 2:K)'];
     Error = norm((X0-X1), 'fro');
     Error = Error/norm(X0, 'fro'); 
     
     err = norm(x1-x0)/norm(x0);
     X0 = X1;
     x0 = x1;
     It = It + 1;    
  end
  
% output
  x = x1;

end

