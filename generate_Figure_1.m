

%% Figure 1 generator
% Script to generate Figure 1 in the submitted paper
% Figure 1.(a) plot the behaviour of F_rho(X^v) by pMAP with fixed rho 
%   while solving the time series denoising problem

%% Setting environment

clear all; clc
AddPath();

fprintf('\n-------------------------------------------------');
fprintf('\n Start generating Figure 1 for pMAP paper');

rng('default') % using default random seed
format short

%% generate the result kpi (functional value and Xerror gap) for pMAP
%   N is the data length
%   r = objective rank * 0.5
%   theta is the noise level of the observed data

N = 2000;
r = 10;
theta = 0.2;

infos_result = result_generator_figure_1(N, r, theta);

%% Generate Figure 1 using the result obtained from result_generator_figure_1()
% Fig.1(a): plot of functional value F_rho at each iteration
% Fig.1(b): plot of Xerror || X^{nu} - X^{nu+1} || at each iteration

figure_generator_figure_1(infos_result)
