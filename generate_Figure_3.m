
%%  Figure 3 generator
%   Script to generate Figure 3 in submitted paper
%   This figure shows the coefficients reconstruction results by FIHT, PGD
%   and pMAP from a random incomplete spectral sparse signal 

%% Setting environment


clear all;clc
AddPath();

fprintf('\n-------------------------------------------------');
fprintf('\n Start generating Figure 3 for pMAP paper');

rng('default')

recovered_series = result_generator_figure_3();

%% Plot Figure 3
%  Figure 1: FIHT, r = 20 (3.a)
%  Figure 2: PGD, r = 20 (3.c)
%  Figure 3: pMAP, r = 20 (3.e)
%  Figure 4: FIHT, r = 40 (3.b)
%  Figure 5: PGD, r = 40 (3.d)
%  Figure 6: pMAP, r = 40 (3.f)

figure_generator_figure_3(recovered_series);


