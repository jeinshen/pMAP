function [pars] = setting_parameters_table_1(weight_choice, N, r)
% Generate necessary parameters for all choices used in Table 1, Figure 1 -
% Figure 2.

L = floor(N/2)+1; % Windows length (numer of rows) of Hankel Matrix
K = N - L + 1; % Number of rows of Hankel Matrix

% Set weight matrix W according to weight choice
% weight_choice == 1: using default SSA-type weight matrix
% weight_choice == 2(or anything else): using equalised weight for all 
%  observations 
if weight_choice == 1
    w = L * ones(N, 1);
    w(1:L) = (1:L);
    w(N-L+1:N) = (L:-1:1);
else
    w = ones(N, 1);
end
W = wmatrix(w, K, L);
pars.W = W/(norm(W, 'fro'));
pars.w = w;

pars.tol = 1.0e-5; % Stop condition: tol
pars.ftol = 1.0e-7; % Stop condition: ftol
pars.gtol = 1.0e-8; % Stop condition: gtol
pars.maxit = 200;   % Stop Condition: max iterative steps
pars.rho = 1e-2/N; % Initialised penalty parameter rho_0
pars.rhoratio = 1.1; % Rho updating ratio
pars.fastswitch = 20000;  % Steps  switch to fast algorithm;

pars.r = r; % Objective rank
pars.L = L;
pars.N = N; % Time series length
pars.K = K;  


end

