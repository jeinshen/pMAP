function [si,iter,ratio,x]=IHT_2D(obs,n1,n2,r,K,maxit,tol,trace)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Iterative hard thresholding algrotihm for 2D spectrally sparse signals.
%
%Inputs
% obs: 1D vector that stores the observed samples. 
% n1: row number of the full signal to be reconstructed.
% n2: column number of the full signal to be reconstructed.
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
% [K,ox,f]=generate_signal([64 32],5,256);obs=ox(K);
% [si,iter,ratio,x]=IHT_2D(obs,64,32,5,K,500,1e-5,1);
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

if mod(n1,2)
    p1=(n1+1)/2;
    d1=[1:p1 p1-1:-1:1]';
else
    p1=n1/2;
    d1=[1:p1 p1 p1-1:-1:1]';
end
if mod(n2,2)
    p2=(n2+1)/2;
    d2=[1:p2 p2-1:-1:1]';
else
    p2=n2/2;
    d2=[1:p2 p2 p2-1:-1:1]';
end

DD=kron(d2,d1);
DD=reshape(DD,n1,n2);

q1=n1+1-p1;
q2=n2+1-p2;

l1=p1*p2;
l2=q1*q2;

epsilon=eps(l1*l2);

ratio=zeros(maxit,1);

% indecies pre-computed for fhmvmultiply_2D to use
ind1=zeros(l2,1);
for i=1:q2
    ind1((i-1)*q1+1:i*q1)=(i-1)*n1+1:(i-1)*n1+q1;
end
ind2=zeros(l1,1);
for i=1:p2
    ind2((i-1)*p1+1:i*p1)=(q2+i-2)*n1+q1:(q2+i-1)*n1;
end
if mod(n1,2) && mod(n2,2)
    ind3=ind1;
    ind4=ind2;
else
    ind3=zeros(l1,1);
    for i=1:p2
        ind3((i-1)*p1+1:i*p1)=(i-1)*n1+1:(i-1)*n1+p1;
    end
    ind4=zeros(l2,1);
    for i=1:q2
        ind4((i-1)*q1+1:i*q1)=(p2+i-2)*n1+p1:(p2+i-1)*n1;
    end
end

% Initialization via one step hard thresholding
N=n1*n2;
m=length(obs); % number of observed samples

alpha=N/m;
b=alpha*obs;

x=zeros(n1,n2); % the full signal to be reconstructed
x(K)=b;

opts=[];opts.eta=1e-16; % set tolerance for lansvd

Yforward=@(z) fhmvmultiply_2D(x,z,q1,q2,ind1,ind2);
Ytranspose=@(z) fhmvmultiply_2D(conj(x),z,p1,p2,ind3,ind4);
[U,S,V]=lansvd(Yforward,Ytranspose,l1,l2,r,'L',opts);
s=diag(S);

x=zeros(n1,n2);
for i=1:r
    ui=reshape(U(:,i),p1,p2);
    vi=reshape(V(:,i),q1,q2);
    x=x+s(i)*conv_fft(ui,conj(vi));
end
x=x./DD;

% Iterations
for iter=1:maxit
    
    x_old=x;
    
    x(K)=b+(1-alpha)*x(K); % gradient descent with step size alpha
    
    Yforward=@(z) fhmvmultiply_2D(x,z,q1,q2,ind1,ind2);
    Ytranspose=@(z) fhmvmultiply_2D(conj(x),z,p1,p2,ind3,ind4);
    [U,S,V]=lansvd(Yforward,Ytranspose,l1,l2,r,'L',opts);
    s=diag(S);
    
    srev=s(r:-1:1);
    r=r-sum(cumsum(srev)/sum(s)<epsilon);
    
    x=zeros(n1,n2);
    for i=1:r
        ui=reshape(U(:,i),p1,p2);
        vi=reshape(V(:,i),q1,q2);
        x=x+s(i)*conv_fft(ui,conj(vi));
    end
    x=x./DD;

    ratio(iter)=norm(x(:)-x_old(:))/norm(x_old(:));
    
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