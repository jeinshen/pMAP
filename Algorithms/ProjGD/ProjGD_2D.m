function [si,iter,x,ratio,fv,gm,step,t] = ... 
    ProjGD_2D(obs,n1,n2,r,K,maxit,tol_1,tol_2,opt_1,opt_2,trace)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projected gradient descent for recovering 2D spectrally sparse signals
% via low rank 2-level Hankel matrix completion.
%
%Inputs
%   obs: vector of observed samples.
%    n1: length of the signal to be recovered in the first dimension.
%    n2: length of the signal to be recovered in the second dimension.
%     r: model order of the spectrally sparse signal.
%     K: locations for the observed samples.
% maxit: maximum number of iterations.
% tol_1: stopping criterion measuring relative change in the signal. 
% tol_2: stopping criterion measuring relative change in function value.
% opt_1: options for projection. 0 (from initialization), >0 (user input).
% opt_2: options for stepsize. 0 (line search), 1 (by minimizing a fourth
%        order polynomail about the stepsize).
% trace: display relative change in signal and function value per iteraiton
%        if not zero.
%
%Outputs
%    si: success indicator, 1 if convergence is achieved, 0 otherwise.
%  iter: iteration number at convergence.
%     x: recovered signal.
% ratio: relative change in signal per iteration.
%    fv: function value per iteration.
%    gm: gradient magnitude per iteration.
%  step: stepsize per iteration. 
%     t: time elapsed per iteration.
%
%References: Cai J F, Wang T, Wei K. Spectral Compressed Sensing via 
%Projected Gradient Descent[J]. arXiv preprint arXiv:1707.09726, 2017.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

si = 0;
ratio = zeros(maxit,1);
fv = zeros(maxit+1,1);
gm = zeros(maxit,1);
step = zeros(maxit,1);
t = zeros(maxit+1,1);

% Initialization via One Step Hard Thresholding
if mod(n1,2)
    p1 = (n1+1)/2;
    d1 = [1:p1 p1-1:-1:1].';
else
    p1 = n1/2;
    d1 = [1:p1 p1 p1-1:-1:1].';
end
if mod(n2,2)
    p2 = (n2+1)/2;
    d2 = [1:p2 p2-1:-1:1].';
else
    p2 = n2/2;
    d2 = [1:p2 p2 p2-1:-1:1].';
end

q1 = n1+1-p1;
q2 = n2+1-p2;

l1 = p1*p2; % number of rows of the 2-level Hankel matrix
l2 = q1*q2; % number of columns of the 2-level Hankel matrix

DD = kron(d2,d1);
DD = reshape(DD,n1,n2);
D = sqrt(DD);

m = numel(obs); % number of observed samples
N = n1*n2;
alpha = N/m;

b = alpha*obs;

x = zeros(n1,n2);
x(K) = b;
x = reshape(x,N,1);

L1 = 2^nextpow2(n1); % next power of 2 for faster fft
L2 = 2^nextpow2(n2); % next power of 2 for faster fft

LL{1} = L1;
LL{2} = L2;

L = LL{1}*LL{2};

Block{1} = 1:n1;
Block{2} = 1:n2;

% indeces for fhmvmultiply to use
ind1 = zeros(l2,1);
for i = 1:q2
    ind1((i-1)*q1+1:i*q1) = (i-1)*n1+1:(i-1)*n1+q1;
end
ind2 = zeros(l1,1);
for i = 1:p2
    ind2((i-1)*p1+1:i*p1) = (q2+i-2)*n1+q1:(q2+i-1)*n1;
end
if mod(n1,2) && mod(n2,2)
    ind3 = ind1;
    ind4 = ind2;
else
    ind3 = zeros(l1,1);
    for i = 1:p2
        ind3((i-1)*p1+1:i*p1) = (i-1)*n1+1:(i-1)*n1+p1;
    end
    ind4 = zeros(l2,1);
    for i = 1:q2
        ind4((i-1)*q1+1:i*q1) = (p2+i-2)*n1+p1:(p2+i-1)*n1;
    end
end

fft_x = fft(x,L);
Yforward = @(z) fhmvmultiply(fft_x,z,N,L,ind1,ind2);
fft_conj_x = fft(conj(x),L);
Ytranspose = @(z) fhmvmultiply(fft_conj_x,z,N,L,ind3,ind4);
opts.eta = 1e-16; % to maintain accuracy
[U,S,V] = lansvd(Yforward,Ytranspose,l1,l2,r,'L',opts);

s = sqrt(diag(S));
s = reshape(s,1,r);

if opt_1 == 0 % projection parameter from initialization
    
    U_2_inf = max(sqrt(sum(U.*conj(U),2)));
    V_2_inf = max(sqrt(sum(V.*conj(V),2)));    
    
    param = max(U_2_inf,V_2_inf)*s(1);    
    
    U = bsxfun(@times,U,s);
    V = bsxfun(@times,V,s);

elseif opt_1 > 0 % user-input projection parameter     
    
    param = opt_1;
    
    U = bsxfun(@times,U,s);
    V = bsxfun(@times,V,s); 
    
    U_r2 = sqrt(sum(U.*conj(U),2));
    V_r2 = sqrt(sum(V.*conj(V),2));   
    
    IND_U = find(U_r2>param);
    IND_V = find(V_r2>param);   
    
    if isempty(IND_U) == 0  
        weight = param./U_r2(IND_U);
        U(IND_U,:) = bsxfun(@times,weight,U(IND_U,:));
    end
    
    if isempty(IND_V) == 0
        weight = param./V_r2(IND_V);
        V(IND_V,:) = bsxfun(@times,weight,V(IND_V,:));    
    end
    
else
    
    error('opt_1: Invalid value!')
    
end

fft_U = zeros(L,r); 
fft_V = zeros(L,r);

x = zeros(n1,n2);
for i = 1:r
    ui = reshape(U(:,i),p1,p2);
    vi = reshape(conj(V(:,i)),q1,q2);
    for dim = 1:2
        ui = fft(ui,LL{dim},dim);
        vi = fft(vi,LL{dim},dim);
    end    
    fft_U(:,i) = reshape(ui,L,1);
    fft_V(:,i) = reshape(vi,L,1);
    ui = ui.*vi;
    for dim = 1:2
        ui = ifft(ui,[],dim);
    end
    ui = ui(Block{:});
    x = x+ui;
end
x = x./DD;
y = D.*x;

UtU = U'*U;
VtV = V'*V;

nobs = D(K).*obs; % weighted observed samples

value = function_evaluation(UtU,VtV,y,alpha,K,nobs);
fv(1) = value;

% other variables to be pre-allocated
HV = zeros(l1,r);
HtU = zeros(l2,r);
fft_Ug = zeros(L,r); 
fft_Vg = zeros(L,r);

t(1) = toc;

% Successive Iterations
for iter = 1:maxit
    
    tic;
    
    x_old = x;

    x(K) = b+(1-alpha)*x(K);
    x = reshape(x,N,1);
    
    fft_x = fft(x,L);
    for i = 1:r
        vi = V(:,i);
        HV(:,i) = fhmvmultiply(fft_x,vi,N,L,ind1,ind2);
    end
        
    S = 1/4*(UtU)+3/4*(VtV);
    Ug = HV-U*S;
    
    fft_conj_x = fft(conj(x),L);
    for i = 1:r
        ui = U(:,i);
        HtU(:,i) = fhmvmultiply(fft_conj_x,ui,N,L,ind3,ind4);
    end
    
    S = 1/4*(VtV)+3/4*(UtU);
    Vg = HtU-V*S;
    
    UgtUg = Ug'*Ug;
    VgtVg = Vg'*Vg;
    
    gm(iter) = sqrt(sum(diag(UgtUg+VgtVg))); % gradient magnitude         
                    
    y1 = zeros(n1,n2);
    for i = 1:r
        vi = reshape(conj(Vg(:,i)),q1,q2);
        ui = reshape(fft_U(:,i),LL{1},LL{2});
        for dim = 1:2
            vi = fft(vi,LL{dim},dim);
        end      
        fft_Vg(:,i) = reshape(vi,L,1);        
        ui = ui.*vi;
        for dim = 1:2
            ui = ifft(ui,[],dim);
        end
        ui = ui(Block{:});
        y1 = y1+ui;
    end
    y1 = y1./D;
    
    y2 = zeros(n1,n2);
    for i = 1:r
        ui = reshape(Ug(:,i),p1,p2);
        for dim = 1:2
            ui = fft(ui,LL{dim},dim);
        end
        vi = reshape(fft_V(:,i),LL{1},LL{2});
        fft_Ug(:,i) = reshape(ui,L,1);      
        ui = ui.*vi;
        for dim = 1:2
            ui = ifft(ui,[],dim);
        end
        ui = ui(Block{:});
        y2 = y2+ui;
    end
    y2 = y2./D;
    
    y3 = zeros(n1,n2);
    for i = 1:r
        ui = reshape(fft_Ug(:,i),LL{1},LL{2});
        vi = reshape(fft_Vg(:,i),LL{1},LL{2});
        ui = ui.*vi;
        for dim = 1:2
            ui = ifft(ui,[],dim);
        end
        ui = ui(Block{:});
        y3 = y3+ui;
    end
    y3 = y3./D;
    
    UtUg = U'*Ug;
    UgtU = UtUg';
    VtVg = V'*Vg;
    VgtV = VtVg';
    
    if opt_2 == 0 % stepsize is chosen via line search   
        
        k = 1;
        
        eta = 1; % initial stepsize for line search
        ls = 0; % number of line search

        while k
           
            y_update = y+(eta*y1+eta*y2+eta^2*y3);
            
            UtU_update = UtU+eta*UtUg+eta*UgtU+eta^2*UgtUg;
            VtV_update = VtV+eta*VtVg+eta*VgtV+eta^2*VgtVg;
            
            value_update = function_evaluation(UtU_update,VtV_update,y_update,alpha,K,nobs);
            
            if value_update < value-1/5*eta*(gm(iter)^2)
                k = 0;
            else
                eta = eta/2;
                ls = ls+1;
            end
            
        end
        
        step(iter) = 1/(2^ls);
        
    elseif opt_2 == 1 % minimizing a fourth order polynomial about stepsize
        
        coeff = zeros(1,4); % coefficients of the derivative of the polynomial
        
        yK = y(K);
        y12 = y1+y2;
        y12K = y12(K);
        y3K = y3(K);
        
        TEMP1 = UgtUg-VgtVg;
        TEMP2 = UgtU+UtUg-VgtV-VtVg;
        TEMP3 = UtU-VtV;
        
        coeff(1) = 2*(sum(sum(conj(UgtUg).*VgtVg))-sum(sum(sum(conj(y3).*y3)))+alpha*(y3K'*y3K))...
            +1/4*(sum(sum(conj(TEMP1).*TEMP1)));
        
        coeff(2) = 3*(sum(sum(conj(UgtUg).*VgtV))+sum(sum(conj(UgtU).*VgtVg))...
            -sum(sum(sum(conj(y12).*y3)))+alpha*y12K'*y3K)+3/8*(sum(sum(conj(TEMP1).*TEMP2)));
        
        coeff(3) = sum(sum(conj(UgtUg).*VtV))+2*sum(sum(conj(UtUg).*VgtV))...
            +sum(sum(conj(UtU).*VgtVg))-sum(sum(sum(conj(y12).*y12)))+alpha*(y12K'*y12K)...
            +2*(sum(sum(conj(UgtU).*VgtV))-sum(sum(sum(conj(y).*y3)))+alpha*(yK-nobs)'*y3K)...
            +1/8*(sum(sum(conj(TEMP2).*TEMP2)))+1/4*(sum(sum(conj(TEMP1).*TEMP3)));
        
        coeff(4) = sum(sum(conj(UgtU).*VtV))+sum(sum(conj(UtU).*VgtV))...
            -sum(sum(sum(conj(y).*y12)))+alpha*(yK-nobs)'*y12K+1/8*(sum(sum(conj(TEMP2).*TEMP3)));
        
        coeff = real(coeff);
        
        ETA = roots(coeff);
        ETA = ETA(imag(ETA)==0); % real stepsize only
        
        num = numel(ETA);      
        if num == 1
            eta = ETA;
        else
            coeff = [coeff./[4 3 2 1] 0];
            [~,j] = min(polyval(coeff,ETA));
            eta = ETA(j);
        end
        
        step(iter) = eta;
        
        y_update = y+(eta*y1+eta*y2+eta^2*y3);
        
        UtU_update = UtU+eta*UtUg+eta*UgtU+eta^2*UgtUg;
        VtV_update = VtV+eta*VtVg+eta*VgtV+eta^2*VgtVg;
        
        value_update = function_evaluation(UtU_update,VtV_update,y_update,alpha,K,nobs);
        
    else
        
        error('opt_2: Invalid value!')
    
    end 
    
    U = U+eta*Ug;
    V = V+eta*Vg;
    
    U_r2 = sqrt(sum(U.*conj(U),2));
    V_r2 = sqrt(sum(V.*conj(V),2));
    
    IND_U = find(U_r2>param);
    IND_V = find(V_r2>param);
    
    if isempty(IND_U) && isempty(IND_V) % no need of projection
        
        UtU = UtU_update;
        VtV = VtV_update;
        
        fft_U = fft_U+eta*fft_Ug;
        fft_V = fft_V+eta*fft_Vg;
        
        y = y_update;
        x = y./D;
        
        value = value_update;
        fv(iter+1) = value;
        
    else
        
        if isempty(IND_U) % no need of projection in U
            fft_U = fft_U+eta*fft_Ug;
            UtU = UtU_update;
        else
            weight = param./U_r2(IND_U);
            U(IND_U,:) = bsxfun(@times,weight,U(IND_U,:));
            for i = 1:r
                ui = reshape(U(:,i),p1,p2);
                for dim = 1:2
                    ui = fft(ui,LL{dim},dim);
                end
                fft_U(:,i) = reshape(ui,L,1);
            end
            UtU = U'*U;
        end
        
        if isempty(IND_V) % no need of projection in V
            fft_V = fft_V+eta*fft_Vg;
            VtV = VtV_update;
        else
            weight = param./V_r2(IND_V);
            V(IND_V,:) = bsxfun(@times,weight,V(IND_V,:));
            for i = 1:r
                vi = reshape(conj(V(:,i)),q1,q2);
                for dim = 1:2
                    vi = fft(vi,LL{dim},dim);
                end
                fft_V(:,i) = reshape(vi,L,1);
            end
            VtV = V'*V;
        end
        
        x = zeros(n1,n2);
        for i = 1:r
            ui = reshape(fft_U(:,i),LL{1},LL{2});
            vi = reshape(fft_V(:,i),LL{1},LL{2});
            ui = ui.*vi;
            for dim = 1:2
                ui = ifft(ui,[],dim);
            end
            ui = ui(Block{:});
            x = x+ui;
        end
        x = x./DD;
        y = D.*x;
        
        value = function_evaluation(UtU,VtV,y,alpha,K,nobs);
        fv(iter+1) = value;
        
    end
    
    ratio(iter) = norm(x(:)-x_old(:))/norm(x_old(:));
    
    t(iter+1) = toc;
    
    if trace
        fprintf('Iteration %4d: relative.change = %.10f, gradient.magnitude = %d \n',iter,ratio(iter),gm(iter))
    end
    
    if ratio(iter) < tol_1 || gm(iter) < tol_2
        si = 1;
        ratio = ratio(1:iter);
        fv = fv(1:iter+1);
        gm = gm(1:iter);
        step = step(1:iter);    
        t = t(1:iter+1);
        return; 
    end
    
end
end % main function ends here

function value = function_evaluation(UtU,VtV,y,alpha,K,nobs)
temp = y(K)-nobs;
value_f = sum(sum(conj(UtU).*VtV))-sum(sum(conj(y).*y))+alpha*(temp'*temp);
TEMP = UtU-VtV;
value_g = 1/8*sum(sum(conj(TEMP).*TEMP));
value = real(value_f+value_g);
end