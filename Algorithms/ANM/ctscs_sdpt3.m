function [x,poles,coeffs,q] = ctscs_sdpt3(xVals,I,n)

T0 = clock;

%%%% ALGORITHM PARAMETERS:
sdpt3_tol = 1e-8; % duality gap tolerance for sdpt3.
prony_tol = 1e-5; % tolerance for estimating the rank when we run Prony
verbose = 0; % set to 1 for printing, 0 for no output
%%%%

m = length(I); % number of measurents
Ic = setdiff(1:n,I); % the indices that are unobserved

%%%% SDPT3 setup:  solves min b'y subject to At'y + C >0 
%%%%
%%%% Our formulation:
%%%% min q(1) subject to [T(q), x; x', q(1)] > 0 and x(i) = xVal(i) for i in I.
%%%%
%%%% To convert complex to real,use 
%%%%           M > 0 iff  [ Re(M), -Im(M); Im(M), Re(M) ] >0
%%%%

if verbose, fprintf('Setting up the SDP Data...\n'); end
    
blk{1,1}='s'; % semidefinite constraint
blk{1,2}=2*n+2; % size of semidefinite constraint

At{1} = sparse(nchoosek(2*n+3,2),4*n-2*m-1);
b = zeros(4*n-2*m-1,1);

% Helper matrices used below:
E = speye(n);
E2 = speye(2);
PX = sparse([0,-1; 1 0]);

% q variable:
At{1}(:,1) = svec(blk,speye(2*n+2)); % q(1), just the identity
b(1) = 1;
for j=2:n,
    At{1}(:,2*j-2) = svec(blk,kron(E2,blkdiag(toeplitz(E(:,j)),0))); % real part
    At{1}(:,2*j-1) = svec(blk,kron(PX,blkdiag(toeplitz(-E(:,j),E(:,j)),0))); % imaginary part
end

% x variable:
for j=1:(n-m)
    At{1}(:,2*n+2*j-2) = svec(blk,kron(E2,sparse([Ic(j);n+1],[n+1;Ic(j)],[1;1],n+1,n+1))); % real part
    At{1}(:,2*n+2*j-1) = svec(blk,kron(PX,sparse([Ic(j);n+1],[n+1;Ic(j)],[1;-1],n+1,n+1))); % imaginary part
end

% C is the matrix of given x values:
CR = sparse([I'; (n+1)*ones(m,1)],[(n+1)*ones(m,1); I'],[real(xVals);real(xVals)],n+1,n+1); % real part
CI = sparse([I'; (n+1)*ones(m,1)],[(n+1)*ones(m,1); I'],[imag(xVals);-imag(xVals)],n+1,n+1); % imaginary party
C = kron(E2,CR)+kron(PX,CI);

% Call SDPT3
if verbose, fprintf('Calling sdpt3...\n'); end

opts = sqlparameters;
if verbose, opts.printlevel = 3; 
else, opts.printlevel = 0;
end
opts.gaptol = sdpt3_tol;

[obj,X,y,Z,info,runhist]=sqlp(blk,At,C,b,opts);

% Extract x:
x = zeros(n,1);
x(I) = xVals;
idxR = 2*n - 2 + 2*(1:(n-m));
idxC = 2*n - 1 + 2*(1:(n-m));
x(Ic) = -y(idxR)-i*y(idxC);

% Extract q:
idxR = [1,2*(2:n)-2];
idxC = 2*(2:n)-1;
q = -y(idxR) -i*[0;y(idxC)];

%%%% Prony's techniqe:
if verbose, fprintf('Running Prony...\n'); end
% Estimate the rank of T(q):
[U,E]=eig(toeplitz(q));
e = diag(E);
order_guess = min(sum(e>prony_tol),n-1);

% With this, set up a companion system, and solve it to find the poles:
X0 = hankel(q(1:n-order_guess),q(n-order_guess:n-1));
X1 = hankel(q(2:n-order_guess+1),q(n-order_guess+1:n));
[UX1,SX1,VX1]=svd(X1,'econ');
s = diag(SX1);
[V,D] = eig(VX1*diag(1./s)*UX1'*X0);

poles = sort(mod(atan2(imag(diag(D)),real(diag(D))),2*pi));

% From the poles, solve a linear system for the coefficients:
U2 = exp( 1i*(0:(n-1))'*poles');
coeffs = U2(I,:)\xVals;

if verbose, fprintf('total time %.2f s\n',etime(clock,T0)); end
