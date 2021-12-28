

%%  Figure 2 generator
%   Script to generate Figure 2 in the submitted paper;
%   This figure shows the behviour of sigma_{r+1}/sigma_r at each iteration
%   for pMAP for time series denoising problems, setting rho is updated at 
%   each iteration in Figure 2(a) and fixed without updating in Figure 2(b)

%% Setting environment

clear all; clc
AddPath();

fprintf('\n-------------------------------------------------');
fprintf('\n Start generating Figure 2 for pMAP paper');

rng('default') % using default random seed
format short

%% generate the result kpi (sigma_{r+1}/sigma_r) with increasing or fixed rho for pMAP
%   N is the data length
%   r = objective rank * 0.5
%   theta is the noise level of the observed data

N = 2000;
r = 10;
theta = 0.2;

infos_result = result_generator_figure_2(N, r, theta);

%% Generate Figure 2: 
%   Fig.2(a): plot of sigma_r+1 / sigma_r at each iteration rho will be 
%    updated (increasing) at each iterate
%   Fig.2(b): plot of functional value F_rho at each iteration rho will be 
%    fixed and not updated 

figure_generator_figure_2(infos_result)
