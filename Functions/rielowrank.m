function [ X, x, U, V ] = rielowrank( x, U, V, r)
%rielowrank This function aims to calculate the low rank projection of a
%           given Hankel matrix, embedded by vector x; This method is
%           inspired by Riemannian optimization.
%
%Reference: Cai J F, Wang T, Wei K. Fast and provable algorithms for 
%           spectrally sparse signal reconstruction via low-rank Hankel 
%           matrix completion[J]. Applied and Computational Harmonic 
%           Analysis, 2019, 46(1): 94-121.
%
%Input      x: estimated signal from last iterate.
%           U: contains the left singular vectors from the SVD in last 
%               iterate. H(x) = U*S*V';
%           V: contains the right singular vectors from the SVD in last 
%               iterate. H(x) = U*S*V';
%           r: objective rank
%
%Output     X: Low rank matrix obtained in current iterate
%           x: estimated signal from this iterate
%           U: updated left singular vectors from the SVD in current 
%               iterate. H(x) = U*S*V';
%           V: updated right singular vectors from the SVD in current 
%               iterate. H(x) = U*S*V';


N = length(x);
L = size(U, 1);
K = N - L + 1;
epsilon = eps(L * K);

if mod(N, 2)  % if N is odd
    
    % Set necessary variables
    
    DD=[1:L L-1:-1:1]';
    HxUc=zeros(L,r);
    Z=zeros(r,r);
    Uc=conj(U);
    
    for i=1:r
        uic=Uc(:,i);
        HxUc(:,i)=fhmvmultiply_1D(x,uic); % To obtain H(x_new)*Uc'
    end
    
    C=(U')*HxUc; 
    Xt=HxUc-U*C; 
    [Q,R]=qr(Xt,0); % qr decomposition
    M=[C R.';R Z];
    
    [Uc,S,Vc]=svd(M);
    s=diag(S(1:r,1:r)); % svd on M;
    
    %srev=s(r:-1:1);
    %r=r-sum(cumsum(srev)/sum(s)<epsilon);
    
    Uc=Uc(:,1:r);
    Vc=Vc(:,1:r);
    
    gamma=-(angle(Uc)+angle(Vc));
    Uc=Uc.*sqrt(exp(1i*gamma));
    U=[U Q]*Uc;
    
    x_new=zeros(N,1);
    for i=1:r
        ui=U(:,i);
        x_new=x_new+s(i)*conv_fft(ui,ui);
    end
    x_new=x_new./DD;

    X = U(:, 1:r) * S(1:r, 1:r) * conj(U(:, 1:r))'; % Low rank matrix; Not Hankel
    x = x_new;
    
else % if N is even
    
    DD=[1:L L L-1:-1:1]';
    UtHx=zeros(r,K);
    HxV=zeros(L,r);
    Z=zeros(r,r);
    for i=1:r
        ui=U(:,i);
        UtHx(i,:)=(fhmvmultiply_1D(conj(x),ui))';
        vi=V(:,i);
        HxV(:,i)=fhmvmultiply_1D(x,vi);
    end
    
    C=UtHx*V;
    Xt=UtHx-C*(V');
    Xc=Xt';
    Yc=HxV-U*C;
    [Q1,R1]=qr(Xc,0);
    [Q2,R2]=qr(Yc,0);
    M=[C R1';R2 Z];
    
    [Uc,S,Vc]=svd(M);
    s=diag(S(1:r,1:r));
    
    srev=s(r:-1:1);
    r=r-sum(cumsum(srev)/sum(s)<epsilon);
    
    U=[U Q2]*Uc(:,1:r);
    V=[V Q1]*Vc(:,1:r);
    
    x_new=zeros(N,1);
    for i=1:r
        ui=U(:,i);
        vi=V(:,i);
        x_new=x_new+s(i)*conv_fft(ui,conj(vi));
    end
    
    x_new=x_new./DD;
    X = U(:, 1:r) * S(1:r, 1:r) * conj(V(:, 1:r))'; % Low rank matrix; Not Hankel
    x = x_new;
end

    
end