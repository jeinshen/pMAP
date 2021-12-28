function kpi_table = experiment_table_2_2_anm(instance_amount)
% Main function to generate and return the kpi table used in the second half
%  of Table 2 of submitted paper, in case of 60% data in a signal are
%  observed
% Note this function is used for ANM solver only

% Create empty kpi table
kpi_table = NaN(8, 1);
solver_names = "ANM";

fprintf('\n-------------------------------------------------');
fprintf('\n Recovering incomplete signal without noises');
fprintf('\n Table 2 when m/n = 0.6, for ANM only');
fprintf('\n-------------------------------------------------');

for data_set = 1 : 2 % Iterate over all data_set
    
    switch data_set
        case 1
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.1, n/m/r = 499/300/20');
            r = 20;
            Y = complex_signal_data_generator(499, 300, r, 0, instance_amount);
        case 2
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.2, n/m/r = 499/300/40');
            r = 40;
            Y = complex_signal_data_generator(499, 300, r, 0, instance_amount);
    end
   
    Kset = Y{1};   % Kset: contains the location of known observations
    OXset = Y{2};  % OXset: contains the original noiseless signal
    
    for solver_choice = 1 : 1
        
        fprintf('\n Solving the signal recovery problem via Solver %s', solver_names(solver_choice))

        % Create kpi table for a specific solver
        kpi_table_solver = zeros(instance_amount, 3);

        for instance = 1 : instance_amount % Iterate over all 50 instances
            K1 = Kset(:, instance);
            ox_signal = OXset(:, instance);
            pars = setting_parameters_table_2(r, K1, ox_signal);
            [iterate, time, error] = solve_problem_table_2_anm(pars);
            kpi_table_solver(instance, :) = [iterate, time, error];
        end
        
        fprintf('\n Solver %s process finished', solver_names(solver_choice));
        
        % Record kpis to main kpi table
        kpi_index = (data_set-1) * 4 + 1;
        % KPI: total number of iterations
        kpi_table(kpi_index, solver_choice) = mean(kpi_table_solver(:,1));
        % KPI: cpu time
        kpi_table(kpi_index + 1, solver_choice) = mean(kpi_table_solver(:,2));
        % KPI: RMSE
        kpi_table(kpi_index + 2, solver_choice) = mean(kpi_table_solver(:,3));
        % KPI: Success recovery rate
        kpi_table(kpi_index + 3, solver_choice) = sum(kpi_table_solver(:,3) < 1.0e-3) / instance_amount;

    end
   
fprintf('\n-------------------------------------------------');

end

end

