function kpi_list = result_generator_figure_4(instance_amount)
% Main function to generate the kpi results used in Figure 4 of sumitted
%   paper.
%   Input data: 3999_1200.mat, including 50 instances of randomly generated
%               imcomplete spectral sparse signal with length 3999 and 1200
%               known observations for each instance.
%   Return: kpi_list, contains the success recover rate of three candidate
%           solvers: FIHT, PGD and pMAP

% create the empty kpi table

fprintf('\n-------------------------------------------------');
fprintf('\n Generate approximation result for Figure 4');

kpi_list = zeros(3, 9);

fprintf('\n-------------------------------------------------');
fprintf('\n Randomly generate test data');

rank = 15;
Y = complex_signal_data_generator(3999, 1200, rank, 0, instance_amount);


Kset = Y{1};
OXset = Y{2};
    
% main loop to iterate all instances for all solvers

fprintf('\n-------------------------------------------------');

for r = 1 : 9
    
    estimate_rank = r * 3 + 3; % estimated rank, range from 3 to 30;
    
    fprintf('\n Solving the signal recovery problem setting rank as %d', estimate_rank)
    
    for solver_choice = 1 : 3
        kpi_table_solver = zeros(50, 1);
        
        for instance = 1 : 50
            K1 = Kset(:, instance); % locations of each known observations
            ox_signal = OXset(:, instance); % original signal
            % the parameter setting in this experiment is the same as Table
            % 2
            pars = setting_parameters_table_2(estimate_rank, K1, ox_signal);
            [~, ~, error] = solve_problem_figure_4(solver_choice, pars);
            kpi_table_solver(instance,1) = error;
        end
        
        kpi_list(solver_choice, r) = sum(kpi_table_solver < 1.0e-3) / 50;
        
    end
    
end

fprintf('\n-------------------------------------------------');
fprintf('\n Result used in Figure 4 are generated');

end

