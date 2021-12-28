%% Script to generate Figure 5 in submitted paper
%  Figure 5 compares the performance of three solvers (FIHT, PGD and pMAP)
%   in incomplete spectral sparse signal recovery when rank is over
%   estimated.


%% Get Data Used For Plotting
clear all; clc
AddPath();

rng('default') % using default random seed
format short

fprintf('\n-------------------------------------------------');
fprintf('\n Start generating Figure 5 for pMAP paper');
over_estimate_rank_kpis = result_generator_figure_5();

%% Plot Figure 5.a (X error plot for three solvers)
figure_generator_figure_5_a(over_estimate_rank_kpis{1});


%% Plot Figure 5.b (final singular value plot for three solvers)
figure_generator_figure_5_b(over_estimate_rank_kpis{2});
