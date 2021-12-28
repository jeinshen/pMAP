function [ x ] = Initial(obs, n, K, r)
%Initial    This is the function for initializing the signal for algorithm 
%           input, assuming that only a subset in a time series is
%           observed.
%    Reference: Cai, J. F., Wang, T., & Wei, K. (2019). Fast and provable 
%               algorithms for spectrallu sparse signal reconstruction via 
%               low rank Hankel matrix completion.
%    Input:
%    Output:



if mod(n,2)
    p=(n+1)/2;
    DD=[1:p p-1:-1:1]';
else
    p=n/2;
    DD=[1:p p p-1:-1:1]';
end

q=n+1-p;

% Initialization via one step hard thresholding
m=length(obs); % number of observed samples

alpha=n/m;
b=alpha*obs;

x=zeros(n,1); % the full signal to be reconstructed
x(K)=b;

opts=[];opts.eta=1e-16; % set tolerance for lansvd

if mod(n,2)

    Yforward=@(z) fhmvmultiply_1D(x,z);
    Ytranspose=@(z) fhmvmultiply_1D(conj(x),z);
    [U,S,V]=lansvd(Yforward,Ytranspose,p,q,r,'L',opts);
    s=diag(S);
    
    gamma=-(angle(U(1,:))+angle(V(1,:)))/2;
    if r==1
        U=U*exp(1i*gamma);
    else
        U=U*sparse(1:r,1:r,exp(1i*gamma));
    end
    
    x=zeros(n,1);
    for i=1:r
        ui=U(:,i);
        x=x+s(i)*conv_fft(ui,ui);
    end
    x=x./DD;
    
    
    else
    
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
end

end

