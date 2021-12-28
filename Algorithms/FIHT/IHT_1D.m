function [si,iter,ratio,x]=IHT_1D(obs,n,r,K,maxit,tol,trace)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Iterative hard thresholding algrotihm for 1D spectrally sparse signals.
%
%Inputs
% obs: 1D vector that stores the observed samples. 
% n: number of the full length of the signal to be reconstructed.
% r: model order of the signal.
% K: indecies for the observed samples in the original signal, obs and K
% should be correspoding to each other.
% maxit: maximum number of iterations allowed (suggested value: 500).
% tol: the iteration scheme is terminated when the relative ratio between 
% successive iterates fall below this value (suggested value: 1e-5).
% trace (optional): if not zero, the relative error after each iteation
% is displayed.
% 
%Outputs
% si: indicator of success, 1 is returned when the relative ratio falls 
% below tol within maxit iterations.
% iter: number of iterations executed.
% ratio: vector that stores relative mean square error after each iteration.
% x: reconstructed signal.
%
%Example
% [K,ox,f]=generate_signal(1024,5,128);obs=ox(K);
% [si,iter,ratio,x]=IHT_1D(obs,1024,5,K,500,1e-5,1);
%
%Reference: Cai, J. F., Wang, T., & Wei, K. (2016). Fast and Provable 
%Algorithms for Spectrally Sparse Signal Reconstruction via Low-Rank Hankel
%Matrix Completion. arXiv preprint arXiv:1606.01567.
%
%Last modified: 15-March-2017
%Please email tianming-wang@uiowa.edu for bug report.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('trace','var')
    trace=0;
end

si=0; % success indicator

if mod(n,2)
    p=(n+1)/2;
    DD=[1:p p-1:-1:1]';
else
    p=n/2;
    DD=[1:p p p-1:-1:1]';
end

q=n+1-p;

epsilon=eps(p*q);

ratio=zeros(maxit,1);

% Initialization via one step hard thresholding
m=length(obs); % number of observed samples

alpha=n/m;
b=alpha*obs;

x=zeros(n,1); % the full signal to be reconstructed
x(K)=b;

opts=[];opts.eta=1e-16; % set tolerance for lansvd

Yforward=@(z) fhmvmultiply_1D(x,z);
Ytranspose=@(z) fhmvmultiply_1D(conj(x),z);
[U,S,V]=lansvd(Yforward,Ytranspose,p,q,r,'L',opts);
s=diag(S);

x=zeros(n,1);
for i=1:r
    ui=U(:,i);
    vi=V(:,i);
    x=x+s(i)*conv_fft(ui,conj(vi));
end
x=x./DD;

% Iterations
for iter=1:maxit
    
    x_old=x;
    
    x(K)=b+(1-alpha)*x(K); % gradient descent with step size alpha
    
    Yforward=@(z) fhmvmultiply_1D(x,z);
    Ytranspose=@(z) fhmvmultiply_1D(conj(x),z);
    [U,S,V]=lansvd(Yforward,Ytranspose,p,q,r,'L',opts);

    s=diag(S);
    
    srev=s(r:-1:1);
    r=r-sum(cumsum(srev)/sum(s)<epsilon);

    x=zeros(n,1);
    for i=1:r
        ui=U(:,i);
        vi=V(:,i);
        x=x+s(i)*conv_fft(ui,conj(vi));
    end
    x=x./DD;
    
    ratio(iter)=norm(x-x_old)/norm(x_old);
    
    if trace
        fprintf('Iteration %4d, RMSE = %.10f \n',iter,ratio(iter))
    end
    
    if ratio(iter)<tol
        si=1;
        ratio=ratio(1:iter);
        break;
    elseif ratio(iter)>1
        break;
    end
    
end