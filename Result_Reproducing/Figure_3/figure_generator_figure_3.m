function figure_generator_figure_3(recovered_series)
% Function to generate Figure 3 in submitted paper

%% Plot Figure from 3.a to 3.f

fprintf('\n-------------------------------------------------');
fprintf('\n Step 3/3: generate and save figures for Figure 3');

for solver_choice = 1 : 3
    
    plot_figure_3_sub_figure_20(solver_choice, recovered_series)
    
    plot_figure_3_sub_figure_40(solver_choice, recovered_series)

end

fprintf('\n-------------------------------------------------');
fprintf('\n Figure 3 succesfully generated. Please find the Fig_3(a).eps - Fig_3(f).eps in /Output folder');
fprintf('\n-------------------------------------------------\n');

end

