function [kpi_table] = experiment_table_1(instance_amount)
% Main function to generate and return the kpi table used in Table 1
%  - input data: this function use cos_simulate_data_geneator() method to
%  randomly generate noiseless and noisy time series data
%  - solver: there are four solvers tested in this experiment, including
%  pMAP, Cadzow and DRI
%  - kpi: there are four main kpis reported in Table 1, including total number
%  of iterations, cpu time, RMSE and ratios of instances closing to best
%  possible RMSEs.
%  For detailed information and references, please see README and submiited
%  paper

% Create empty kpi table
kpi_table = NaN(48, 3);
solver_names = ["pMAP", "Cadzow", "DRI"];
% threshold of close RMSE, used to check if a obtained approximation is
% close enough to the best available approximation.
% For example, 0.05 means if the approximated RMSE is within 5% lower or
% higher than the best RMSE, we say they are 'close'.
threshold_close_RMSE = 0.10;

fprintf('\n-------------------------------------------------');
fprintf('\n Recovering incomplete signal without noises');
fprintf('\n Table 1');
fprintf('\n Results are the average of %d instances', instance_amount);
fprintf('\n-------------------------------------------------');

for data_set = 1 : 1 % Iterate over all data_set
    
    % Randomly generate testing data set
    switch data_set
        case 1
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.1, n/r/theta = 1000/10/0.1');
            N = 1000;       % length of time series
            r = 10;         % objective rank
            theta = 0.1;    % noise level
            [Yn, y] = cos_simulate_data_geneator(N, r/2, theta, instance_amount);
        case 2
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.2, n/r/theta = 1000/10/0.2');
            N = 1000;
            r = 10;
            theta = 0.2;
            [Yn, y] = cos_simulate_data_geneator(N, r/2, theta, instance_amount);
        case 3
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.3, n/r/theta = 1000/10/0.5');
            N = 1000;
            r = 10;
            theta = 0.5;
            [Yn, y] = cos_simulate_data_geneator(N, r/2, theta, instance_amount);
        case 4
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.4, n/r/theta = 2000/20/0.1');
            N = 2000;
            r = 20;
            theta = 0.1;
            [Yn, y] = cos_simulate_data_geneator(N, r/2, theta, instance_amount);
        case 5
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.5, n/r/theta = 2000/20/0.2');
            N = 2000;
            r = 20;
            theta = 0.2;
            [Yn, y] = cos_simulate_data_geneator(N, r/2, theta, instance_amount);
        case 6
            fprintf('\n-------------------------------------------------');
            fprintf('\n Solving dataset No.6, n/r/theta = 2000/20/0.5');
            N = 2000;
            r = 20;
            theta = 0.5;
            [Yn, y] = cos_simulate_data_geneator(N, r/2, theta, instance_amount);
    end
    
    % We consider two different weight choices
    % weight_choice == 1 indicates we use SSA type weight.
    % weight_choice == 2 indicates we use equal weights.
    % details of weight information can be seen from Section 5.1.1 of submitted paper
    
    for weight_choice = 1 : 2 
        
        pMAP_RMSE_achieved = zeros(instance_amount, 1);
        
        pars = setting_parameters_table_1(weight_choice, N, r);
            
        for solver_choice = 1 : 3 % Three candidate solvers are considered
            
            fprintf('\n Solving the signal recovery problem via Solver %s, weight choice is %d', solver_names(solver_choice), weight_choice)

            % 2nd weight choice for solver DRI is skipped because DRI does
            %   not allow flexible weight choice.
            if weight_choice == 2 && solver_choice == 2
                continue
            end
            
            kpi_table_solver = zeros(instance_amount, 3);
            RMSE_gap_to_pMAP = zeros(instance_amount, 1);
            
            for instance = 1 : instance_amount % Iterate over all instances
                yn = Yn(:, instance);
                [iterate, time, error] = solve_problem_table_1(yn, y, ...
                    solver_choice, pars);
                kpi_table_solver(instance, :) = [iterate, time, error];
                % we record RMSE of pMAP for calculating the ratio of
                % instances closing to the best RMSE (pMAP RMSE) for each
                % solver
                if solver_choice == 1
                    pMAP_RMSE_achieved(instance, 1) = error;
                    RMSE_gap_to_pMAP(instance, 1) = 1;
                else
                    RMSE_gap = (error - pMAP_RMSE_achieved(instance, 1)) / pMAP_RMSE_achieved(instance, 1);
                    if RMSE_gap <= threshold_close_RMSE
                        RMSE_gap_to_pMAP(instance, 1) = 1;
                    end
                end
            end
            
            fprintf('\n Solving process for Solver %s and weight choice is %d is finished', solver_names(solver_choice), weight_choice)

            % Record kpis to main kpi table
            kpi_index = (data_set-1)*8 + (weight_choice-1)*4 + 1;
            % KPI: total amount of iterations
            kpi_table(kpi_index, solver_choice) = mean(kpi_table_solver(:,1));
            % KPI: CPU time
            kpi_table(kpi_index + 1, solver_choice) = mean(kpi_table_solver(:,2));
            % KPI: RMSE
            kpi_table(kpi_index + 2, solver_choice) = mean(kpi_table_solver(:,3));
            % KPI: ratio of close approximations
            kpi_table(kpi_index + 3, solver_choice) = sum(RMSE_gap_to_pMAP) / instance_amount;

        end

    end
    
end


end

