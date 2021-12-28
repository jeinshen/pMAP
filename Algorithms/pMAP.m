function [X, pfhist, infos] = pMAP(Y, W, X0, pars, trace)
%%
% penalised Method of Alternating Projection for Time Series of Finite Rank
% by Jian Shen, Jein-Shan Chen, Houduo Qi and Naihua Xiu.
%
% Solve: min f(X) = 0.5 \| W .* (Y-X) \|^2, s.t. X \in H and X \in M_r
% where  W -- weight matrix (LxK)
%        Y -- given Hankel matrix (LxK)
%        H -- subspace of Hankel matrices
%        M_r -- rank r matrix of LxK
%% Input
%
%    W, Y, r as defined in the problem
%    rho -- penalty parameter
%    X0  -- initial LxK matrix (can be taken to be Y)
%
%    pars -- other parameters
%    pars.maxit     -- maximum number of iterations (200)
%    pars.ftol      -- tolerance used for stopping criterion in f (1.0e-3)
%    pars.gtol      -- tolerance used for stopping criterion in g (1.0e-3)
%    pars.rhoratio  -- rho is increased by the ratio
%    pars.rhomode   -- rho is increased evry rhomode steps
%                      (Default: No updating)
%
%   Output
%
%    X      -- final solution
%    pfhist -- history of the penalty objective functional values
%
%    infos  -- onther output
%     infos.fprog  -- the progess in f at the final iterate
%     infos.gerror -- error in g function (g should be 0)
%     infos.iter   -- number of iteration taken
%     infos.fhist  -- history of objective functional values of f
%
%%
if ~isfield(pars, 'maxit')
    maxit = 200;
else
    maxit = pars.maxit;
end

if isfield(pars, 'ftol')
    ftol = pars.ftol;
else
    ftol = 1.0e-4;
end

if isfield(pars, 'gtol')
    gtol = pars.gtol;
else
    gtol = 1.0e-6;
end

if isfield(pars, 'tol')
    tol = pars.tol;
else
    tol = 1.0e-6;
end

if isfield(pars, 'rhotol')
    rhotol = pars.rhotol;
else
    rhotol = 1e3 * min(min(W));
    %rhotol = 1e3 * mean(pars.w);
end

if isfield(pars, 'rhoratio')
    rhoratio = pars.rhoratio;
else
    rhoratio = 1.5;
end

if isfield(pars, 'rhomode')
    rhomode = pars.rhomode;
else
    rhomode = 2*maxit;  % no updating rho
end

if isfield(pars, 'fastswitch')
    fastswitch = pars.fastswitch;
else
    fastswitch = 300;  % no switch to fast algorithm
end

if isfield(pars, 'xerror_plot')
    xerror_plot = pars.xerror_plot;
else
    xerror_plot = 0;  % no xerror_plot computing by default
end

rho = pars.rho;
r = pars.r;
% SVD for the initial point
[U, S, V] = svds(X0, r);
X0r = U*S*V';
X = X0;

% Record the initial objective value
fhist = zeros(maxit, 1);
ghist = fhist;
pfhist = fhist;
Xhist = fhist;
sigmahist = fhist;
xhist = fhist;

f0 = norm( W .* (X0 - Y) , 'fro');
f0 = 0.5*f0^2;
fhist(1)  = f0;
if rho <= rhotol
    rho = rhoratio * rho;
end

% Compute the initial gerror and set fprog
gerror = norm( (X0r - X0), 'fro');
gerror = 0.5*gerror.^2;
pf0 = f0 + rho*gerror;
pfhist(1) = pf0;

% Scale gerror
gerror = (2*gerror)/(norm(X0, 'fro').^2);

% set fprog
if gerror < gtol
    fprog = 0;
else
    fprog = 1;      % to start the algorithm
end
it = 1;

Wrho2 = W.*W + rho;  % W_rho^2
[L, K] = size(W);
idx = hankel(1:L, L:(K-1)+L);  % to be used for computing Hankel Projection

% start the main algorithm

if trace
    fprintf('\n-------------------------------------------------');
    fprintf('\n penalised Method of Alternating Projection');
    fprintf('\n-------------------------------------------------');

    fprintf('\n  it   | rho       fobj       pfobj       fprog      gerror| \n');
    fprintf('[it=%1.0f  | %3.2e  %3.2e   %3.2e    %3.2e   %3.2e|]', ...
        it, rho,  f0,  pf0,  fprog,  gerror);
end
Xerror = 1;


while (it <= maxit) && ( (fprog > ftol) ) && ((gerror > gtol) )  && (Xerror > tol)
    
    A = W.*W.*Y + rho * X0r;
    
    % record x_old from last iteration
    X_old = X;
    
    % The block below is to compute the new iterate X
    sumH = accumarray(idx(:), A(:));
    sumW = accumarray(idx(:), Wrho2(:));
    sumH = sumH ./ sumW;
    X = hankel(sumH(1:L), sumH(L:end));
    
    it = it + 1;
    Xerror1 = norm(X_old-X, 'fro');
    Xerror = Xerror1/norm(X_old, 'fro');
    if xerror_plot == 1
        xerror = norm(W.*(X_old-X), 'fro')/norm(W.*X_old, 'fro');
        xhist(it-1, 1) = xerror;
    end
    
    % Compute the objective function at new X
    f0 = norm( W .* (X - Y) , 'fro');
    f0 = 0.5*f0^2;
    fhist(it+1)  = f0;
    
    % Compute the g function at new X
    if it > fastswitch
        [ X0r, ~, U, V ] = rielowrank(x_new, U, V, r);
    else
        [U, S, V] = svds(X, r+1);
        X0r = U(:, 1:r)*S(1:r, 1:r)*(V(:, 1:r))';
        sigmahist(it, 1) = S(r+1,r+1)/S(r,r);
    end

    X0  = X;
    
    gerror = norm( (X0r - X0), 'fro');
    gerror = 0.5*gerror.^2;
    gerror0 = gerror;        % to be used when rho is increased
    
    % Compute the penalty objective function at X
    pf = f0 + rho*gerror;
    pfhist(it-1) = pf;
    gerror1 = gerror;
    % Scale gerror
    gerror = (2*gerror)/(norm(X0, 'fro').^2);
    
    % Compute fprog
    if pf > pf0
        if trace
            fprintf('\n ');
            disp('the penalty objective is not decreasing');
        end
        fprog = abs(pf0 - pf)/(1+pf0);
        %fprog = 1;
    else
        fprog = abs(pf0 - pf)/(1+pf0);
    end
    pf0 = pf;
    
    
    % Increase rho
    %if mod(it, rhomode) == 0
    if rho <= rhotol
        rho = rhoratio * rho;
        rho = min(rho, 1.0e5);
        Wrho2 = W.*W + rho;
        pf0 = f0 + rho*gerror0;  % adjust pf0 by the new rho
    end
    
    % print out information
    if trace
        fprintf('\n');
        fprintf('[it=%1.0f  | %3.2e  %3.2e   %3.2e    %3.2e   %3.2e|]', ...
            it,  rho,  f0,  pf0,  fprog,  gerror);
    end
    ghist(it-1, 1) = gerror;
    Xhist(it-1, 1) = Xerror1;
end
if trace
    fprintf('\n ')
end

% Output information
infos.pfhist = pfhist;
infos.fhist = fhist;  % sequence of original objective values
infos.fprog = fprog;
infos.it = it;
infos.gerror = ghist;
infos.Xhist = Xhist;
infos.sigmahist = sigmahist;
infos.xhist = xhist;
X = X0;
end

