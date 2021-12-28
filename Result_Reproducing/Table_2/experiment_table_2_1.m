function kpi_table = experiment_table_2_1(instance_amount)
% Main function to generate and return the kpi table used in the first half
%  of Table 2 of submitted paper, in case of 30% data in a signal are
%  observed
%  - input data: this function use complex_signal_data_generator() method to
%  randomly generate incomplete spectral sparese signals
%  - solver: there are five solvers tested in this experiment, including
%  Cadzow, DRI, FIHT, PGD and pMAP
%  - kpi: there are four kpis reported in Table 2, including total number
%  of iterations, cpu time, RMSE and ratios of instances successfully 
%  recovered.
%  For detailed information and references, please see README and submiited
%  paper

% Create empty kpi table
kpi_table = NaN(24, 5);
solver_names = ["Cadzow", "DRI", "FIHT", "PGD", "pMAP"];

fprintf('\n-------------------------------------------------');
fprintf('\n Recovering incomplete signal without noises');
fprintf('\n Table 2 when m/n = 0.3');
fprintf('\n Results are the average of %d instances', instance_amount);
fprintf('\n-------------------------------------------------');


for data_set = 1 : 6 % Iterate over all data_set
    
    % Randomly generate testing data set
    switch data_set
        case 1
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.1, n/m/r = 499/150/10');
            r = 10;
            Y = complex_signal_data_generator(499, 150, r, 0, instance_amount);
        case 2
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.2, n/m/r = 499/150/20');
            r = 20;
            Y = complex_signal_data_generator(499, 150, r, 0, instance_amount);
        case 3
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.3, n/m/r = 999/300/30');
            r = 30;
            Y = complex_signal_data_generator(999, 300, r, 0, instance_amount);
        case 4
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.4, n/m/r = 999/300/40');
            r = 40;
            Y = complex_signal_data_generator(999, 300, r, 0, instance_amount);
        case 5
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.5, n/m/r = 1999/600/60');
            r = 60;
            Y = complex_signal_data_generator(1999, 600, r, 0, instance_amount);
        case 6
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.6, n/m/r = 1999/600/80');
            r = 80;
            Y = complex_signal_data_generator(1999, 600, r, 0, instance_amount);
    end
    
    Kset = Y{1};   % Kset: contains the location of known observations
    OXset = Y{2};  % OXset: contains the original noiseless signal
    
    for solver_choice = 1 : 5
        
        fprintf('\n Solving the signal recovery problem via Solver %s', solver_names(solver_choice))
                
        % Create kpi table for a specific solver
        kpi_table_solver = zeros(instance_amount, 3);

        for instance = 1 : instance_amount % Iterate over all instances
            
            K1 = Kset(:, instance);
            ox_signal = OXset(:, instance);
            pars = setting_parameters_table_2(r, K1, ox_signal);
            [iterate, time, error] = solve_problem_table_2(solver_choice, pars);
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

