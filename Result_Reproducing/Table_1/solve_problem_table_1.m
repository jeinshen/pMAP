function [iterate, time, error] = solve_problem_table_1(yn, y, solver_choice, pars)
% Method to tackle the time series denoising according to the solver choice 
%   and parameters, return main kpis including total numer of
%   iterations (iterate), computing cpu time (time in seconds) and root
%   mean square errors (error).


tStart = cputime;

switch solver_choice
    
    case 1 % pMAP
        L = pars.L;
        A = hankel(yn(1 : L), yn(L : end));
        [X, ~, infos] = pMAP(A, pars.W, A, pars, 0);
        x = [X(1:L, 1); X(L, 2:end)'];
        iterate = infos.it;

    case 2 % Cadzow
        [x, iterate, ~] = Cadzow(yn, pars);

    case 3 % DRI
        [x, iterate, ~, ~] = DRI(yn, pars);        
        
end

time = cputime - tStart;
error = norm(x-y)/sqrt(pars.N);

end

