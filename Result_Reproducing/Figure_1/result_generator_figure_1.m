function [infos_list] = result_generator_figure_1(N, r, theta)
%% fuction to generate the result kpi (infos) for the plotting of Figure 1 
% Input: N is the data length
%        r = objective rank * 0.5
%        theta is the noise level of the observed data
% Output: the infos_list in case W = W_1 and W = W_2, which includes the implementation
%         information of pMAP algorithm.
%% generate the randomly simulated data for figure 1
%    instance_amount is the number of instances we generated for the figure,
%    which is 1 in this experiment.

fprintf('\n-------------------------------------------------');
fprintf('\n Step 1/4: generate random noisy time series data');
instance_amount = 1;

[Yn, ~] = cos_simulate_data_geneator(N, r, theta, instance_amount);
yn = Yn(:, 1);

%% generate kpi when W = W_1
% 1. pars includes all the parameters used for Figure 1.
% 2. weight_choice = 1 means we use default SSA-type weight matrix. Please 
%       see Section 5.1.1 of pMAP paper for detailed information
% 3. pars.rhoratio = 1 means we fix the penalty rho in figure 1, it is
%       because we want to demostrate the functional values
%       and algorithm hebaviours under the same objective function.

fprintf('\n-------------------------------------------------');
fprintf('\n Step 2/4: solve time series denoising problem setting W = W1');
weight_choice = 1;
pars = setting_parameters_table_1(weight_choice, N, r * 2);
pars.rhoratio = 1;

L = pars.L;
A = hankel(yn(1 : L), yn(L : end));
[~, ~, infos] = pMAP(A, pars.W, A, pars, 0);

infos_list{1} = infos;

%% generate kpi when W = W_2
% 1. pars includes all the parameters used for figure 1.
% 2. weight_choice = 2 means we use equalised weight for all observations.
%       Please see Section 5.1.1 of pMAP paper for detailed information
% 3. pars.rhoratio = 1 means we fix the penalty rho in figure 1, it is
%       because we want to demostrate the functional values
%       and algorithm hebaviours under the same objective function.
fprintf('\n-------------------------------------------------');
fprintf('\n Step 3/4: solve time series denoising problem setting W = W2');
weight_choice = 2;
pars = setting_parameters_table_1(weight_choice, N, r * 2);
pars.rhoratio = 1;

[~, ~, infos] = pMAP(A, pars.W, A, pars, 0);

infos_list{2} = infos;

end

