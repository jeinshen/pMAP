function pars = setting_parameters_table_2(r, K1, ox_signal)
% Generate necessary parameters for all choices used in Table 2, Figure 3
% and Figure 4.

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
pars.ox_signal = ox_signal; % Orignal full signal
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
% w_i = 1 if ox_signal is observed, 0 otherwise.
% W is computed via wmatrix method from w
w = zeros(pars.N, 1);
for j = 1 : length(K1)
    w(K1(j, 1)) = 1;
end
W = wmatrix(w, pars.N - pars.L + 1, pars.L);
pars.W = W/(norm(W, 'fro'));
pars.w = w;

end

