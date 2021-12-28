function [iterate, time, error] = solve_problem_figure_4(solver_choice, pars)
% Method to tackle the incomplete signal recovery problem according to the
%   solver choice and parameters, return main kpis including total numer of
%   iterations (iterate), computing cpu time (time in seconds) and root
%   mean square errors (error).

switch solver_choice
        
	case 1 % FHIT
        [~, iterate, ~, x]=FIHT_1D(pars.obs1, pars.N, pars.r, ...
            pars.K1, pars.maxit, pars.tol, 0);
        
	case 2 % PGD
        [~, iterate,x, ~, ~, ~, ~, ~] = ProjGD_1D(pars.obs1, pars.N, ...
            pars.r, pars.K1, pars.maxit,pars.tol,pars.ftol,0,0,0);
        
    case 3 % pMAP
        [X, ~, infos] = pMAP(pars.A, pars.W, pars.A, pars, 0);
        x = [X(1: pars.L, 1); transpose(X(pars.L, 2:end))];
        iterate = infos.it;
end

time = 0;

error = norm(pars.inversObsl.*(x- pars.ox_signal))/norm(pars.inversObsl.* pars.ox_signal);

end

