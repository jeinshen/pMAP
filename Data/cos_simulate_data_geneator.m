

function [Yn, y] = cos_simulate_data_geneator(N, r, theta, instance_amount)
%% Script to generate examples for Experiment 1 (Time series denoising)
% Input: N is the data length
%        r = objective rank * 0.5
%        theta is the noise level of the observed data
%        instance_amount is the number of instances we generated for the figure.
% Output: Yn contains all the instances of randomly generated data
%         y is the noiseless signal series.

%% Set parameters randomly
% The detailed information of d_s, alpha_s, beta_s and tau_s can be found
% fomr the first display in section 5.1.1 of paper: 
%
% A Penalized Method of Alternating Projections for Weighted Low-Rank Hankel 
% Matrix Optimization

ds = rand(r, 1) * 1000;                     % d_s
alphas = 0.999 + 0.002 * rand(r ,1);        % alpha_s
betas =  2 * pi ./ (6 + 12 * rand(r, 1));   % beta_s
taus = (rand(r, 1) - 0.5 ) * 2 * pi;        % tau_s

%% Generate noiseless time series with length N
y = zeros(N, 1);
for i = 1 : N
    y(i,1) = sum(ds .* (alphas .^ i) .* cos( betas .* i - taus));
end

%% Generate noise observations yn = y + n and store in Yn
Yn = zeros(N, instance_amount);
for i = 1 : instance_amount
    n = normrnd(0, 1, N, 1); 
    n = theta * norm(y, 'fro') * n / norm(n, 'fro'); % generate noises and add to noiseless time series
    Yn(:,i ) = y + n;
end

end


