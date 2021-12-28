function recovered_signal = result_generator_figure_3()
% Function to generate the coefficient reconstruction results y FIHT, PGD
%  and pMAP from a random incomplete spectral sparse signal 

% Create empty matrix for recovered signal
recovered_signal = zeros(499, 8);

%% Generate recovered signals by FIHT, PGD and pMAP in case of r = 20

fprintf('\n-------------------------------------------------');
fprintf('\n Step 1/3: Generate recovered signals by FIHT, PGD and pMAP in case of r = 20');

load 499_150_20.mat

Kset = Y{1};        
OXset = Y{2};

K1 = Kset(:, 1);
ox_signal = OXset(:, 1);
pars = setting_parameters_table_2(r, K1, ox_signal);

recovered_signal(:, 1) = ox_signal;

[~, ~, ~, x_fiht]=FIHT_1D(pars.obs1, pars.N, pars.r, ... 
    pars.K1, pars.maxit, pars.tol, 0);
recovered_signal(:, 2) = x_fiht;

[~, ~,x_pgd, ~, ~, ~, ~, ~] = ProjGD_1D(pars.obs1, pars.N, ...
    pars.r, pars.K1, pars.maxit,pars.tol,pars.ftol,0,0,0);
recovered_signal(:, 3) = x_pgd;

[X, ~, ~] = pMAP(pars.A, pars.W, pars.A, pars, 0);
x_pmap = [X(1: pars.L, 1); transpose(X(pars.L, 2:end))];
recovered_signal(:, 4) = x_pmap;

%% Generate recovered signals by FIHT, PGD and pMAP in case of r = 40

fprintf('\n-------------------------------------------------');
fprintf('\n Step 2/3: Generate recovered signals by FIHT, PGD and pMAP in case of r = 40');

load 499_300_40.mat
 
Kset = Y{1};        % Input Data
OXset = Y{2};
 
K1 = Kset(:, 1);
ox_signal = OXset(:, 1);
pars = setting_parameters_table_2(r, K1, ox_signal);
 
recovered_signal(:, 5) = ox_signal;
 
[~, ~, ~, x_fiht]=FIHT_1D(pars.obs1, pars.N, pars.r, ... 
    pars.K1, pars.maxit, pars.tol, 0);
recovered_signal(:, 6) = x_fiht;

[~, ~,x_pgd, ~, ~, ~, ~, ~] = ProjGD_1D(pars.obs1, pars.N, ...
    pars.r, pars.K1, pars.maxit,pars.tol,pars.ftol,0,0,0);
recovered_signal(:, 7) = x_pgd;
 
[X, ~, ~] = pMAP(pars.A, pars.W, pars.A, pars, 0);
x_pmap = [X(1: pars.L, 1); transpose(X(pars.L, 2:end))];
recovered_signal(:, 8) = x_pmap;


end
