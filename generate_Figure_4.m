
%% Script to generate Figure 4 in sumitted paper
%  This figure shows the performance of candidate solvers in spectral
%   sparse signal recovering when the input rank is misestimated.
%  Candidate solver including FIHT, PGD and pMAP.
%  For detailed information and references, please see README


%% Generate the approximation results used in Figure 4

clear all; clc
AddPath();
rng('default') % using default random seed
format short

instance_amount = 50;

fprintf('\n-------------------------------------------------');
fprintf('\n Start generating Figure 4 for pMAP paper');

kpi_data_incorrect_rank = result_generator_figure_4(instance_amount);

%% Generate Figure 4 using the generated results

figure_generator_figure_4(kpi_data_incorrect_rank)
