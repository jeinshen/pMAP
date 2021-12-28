function [x,poles,coeffs,x2] = ctscs_dast(xVals,I,n,upsample)

%ctscs_dast(xVals,I,n,upsample)
T0 = clock;
%%%% ALGORITHM PARAMETERS:
coeff_tol = 1e-3; % tolerance where we clip the support
verbose = 0; % set to 1 for printing, 0 for no output
%%%%

N = 2^ceil(log2(upsample*n));
m = length(I); % number of measurents

F = N*ifft(eye(N));
FI = F(I,:);
cvx_precision best;
cvx_solver sdpt3;
cvx_begin quiet
variable c(N) complex;
minimize norm(c,1)
subject to
xVals == FI*c;
cvx_end

x = F(1:n,:)*c;
idx = find(abs(c)>coeff_tol*max(abs(c)));
poles = 2*pi*(idx-1)/N; % poles are on a grid as we used the fft of size N
coeffs = F(I,idx)\xVals;
x2 = F(1:n,idx)*coeffs;

if verbose, fprintf('total time %.2f s\n',etime(clock,T0)); end