function pars = setting_parameters_table_3(r, K1, ox_signal, ox1)
% Generate necessary parameters for all choices used in Table 3

pars.N = length(ox_signal); % Signal length
pars.L = floor((pars.N+1)/2); % Windows length (numer of rows) of Hankel Matrix
pars.r = r; % Objective rank

pars.fastswitch = 2000; % Steps switch to fast algorithm;
pars.tol = 1.0e-5; % Stop condition: tol
pars.ftol = 1.0e-7; % Stop condition: ftol
pars.gtol = 1.0e-8; % Stop condition: gtol
pars.maxit = 200;   % Stop Condition: max iterative steps

% Initialised penalty parameter rho_0
pars.rho = (1.0e-2 * length(K1) )/(pars.N * pars.N);

pars.K1 = K1; % Locations (index) of known observations
pars.ox1 = ox1; % Orignal full signal
pars.ox_signal = ox_signal; % Signal with noisey observations
pars.obs1=ox_signal(K1); % Incomplete observations
pars.M = length(K1); % Number of known observations

% Vector to store the location of known observations.
% 1 is the data is known, 0 otherwise
ObsL = zeros(pars.N, 1);
for j = 1 : length(K1)
	ObsL(K1(j, 1)) = 1;
end
pars.ObsL = ObsL;

% Vector to store the location of unknown data.
% 0 is the data is known, 1 otherwise
pars.inversObsl = ones(pars.N, 1) - ObsL;

% Construct initial input Hankel matrix
x1 = Initial(pars.obs1, pars.N, pars.K1, pars.r);
pars.Yn = x1 .* pars.inversObsl + ox_signal .* pars.ObsL;
pars.A = hankel(pars.Yn(1:pars.L), pars.Yn(pars.L:end));

% Set weight matrix W from the incomplete signal.
% w_i = 1 if ox_signal is observed and noisy, 100 if observation is 
%  noiseless, 0 if the data is unobserved.
% W is computed via wmatrix method from w
w = zeros(pars.N, 1);
for j = 1 : floor(length(K1) *1/3)
    w(K1(j, 1)) = 1;
end
for j = (floor(length(K1) *1/3)+1) : length(K1)
    w(K1(j, 1)) = 100;
end

W = wmatrix(w, pars.N - pars.L + 1, pars.L);
pars.W = W/(norm(W, 'fro'));
pars.w = w;

end

