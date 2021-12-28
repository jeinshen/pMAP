clear;clc;
n = 30;
k = 4;
nn = 0:n-1;
nn = nn(:);
f0 = rand(k,1);
A = exp(1i*2*pi*nn*f0');
c0 = rand(k,1);%randn(k,1) + 1i*randn(k,1);

x = A*c0;

%%
cvx_begin quiet
variable u0;
variable u(n-1) complex;
variable Tp(n,n) hermitian;
variable t;
minimize ((u0 + t)/2)
subject to
Tp == toeplitz([u0;conj(u)],[u0;conj(u)]');
[Tp x;
    x' t] == semidefinite(n+1,n+1);
cvx_end

fprintf('cvx_status = %s, opt = %8.6f\n', cvx_status, (u0+t)/2);

%%
N = 4000;
ff = (0:N-1)/N;
ff = ff(:);


AN = exp(1i*2*pi*nn*ff');

cvx_begin quiet
variable ch(N) complex;
minimize norm(ch,1)
subject to
x == AN*ch;
cvx_end

fprintf('cvx_status = %s, opt = %8.6f\n', cvx_status, norm(ch,1));

%%
N = 2*n+1;
ff = sort(rand(N,1));
ff = ff(:);
A0 = exp(1i*2*pi*nn*ff');
cvx_begin quiet
variable p(N);
variable t;
minimize ((t + sum(p))/2)
subject to
%p >= 0;
[A0*diag(p)*A0' x;
    x' t] == semidefinite(n+1,n+1);
cvx_end

fprintf('cvx_status = %s, opt = %8.6f\n', cvx_status, (t + sum(p))/2);

