function over_estimate_rank_kpis = result_generator_figure_5()
% Function to generate kpi information for plotting figure 5
% Output: over_estimate_rank_kpis{1}: the x error at each iteration
%         over_estimate_rank_kpis{2}: the singular values of final output

fprintf('\n-------------------------------------------------');
fprintf('\n Step 1/4: generate random testing data');

rank = 15;
Y = complex_signal_data_generator(3999, 1200, rank, 0, 50);

Kset = Y{1};
OXset = Y{2};

estimate_rank = 21;

fprintf('\n-------------------------------------------------');
fprintf('\n Step 2/4: solve incomplete signal recovery problem with incorrect rank');

K1 = Kset(:, 2); 
ox_signal = OXset(:, 2); 
pars = setting_parameters_table_2(estimate_rank, K1, ox_signal);
pars.xerror_plot = 1;

[~, ~, ratio_FIHT, x_FIHT]=FIHT_1D(pars.obs1, pars.N, pars.r, pars.K1,...
    pars.maxit, pars.tol, 0);
X_FIHT = hankel(x_FIHT(1:pars.L), x_FIHT(pars.L:end));
singular_values{1} = generate_singular_values(X_FIHT);

[~, ~, x_PGD, ratio_PGD, ~, ~, ~, ~] = ProjGD_1D(pars.obs1, pars.N,pars.r,...
    pars.K1, pars.maxit,pars.tol,pars.ftol,0,0,0);
X_PGD = hankel(x_PGD(1:pars.L), x_PGD(pars.L:end));
singular_values{2} = generate_singular_values(X_PGD);

[X_pMAP, ~, infos] = pMAP(pars.A, pars.W, pars.A, pars, 0);
singular_values{3} = generate_singular_values(X_pMAP);

xerror_ratios{1} = ratio_FIHT;
xerror_ratios{2} = ratio_PGD;
xerror_ratios{3} = infos.xhist;

over_estimate_rank_kpis{1} = xerror_ratios;
over_estimate_rank_kpis{2} = singular_values;
end

function singular_values = generate_singular_values(X)
    [~, S, ~] = svds(X, 50);
    singular_values = diag(S);
end