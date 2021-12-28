function [iterate, time, error] = solve_problem_table_2_anm(pars)
% Method to tackle the incomplete signal recovery problem by ANM according 
%   to parameters, return main kpis including total numer of computing cpu 
%   time (time in seconds) and root mean square errors (error).

tStart = cputime;

[x, ~, ~] = ANM(pars.obs1, pars.K1, pars.N, pars.ox_signal);

time = cputime - tStart;

error = norm(pars.inversObsl.*(x- pars.ox_signal))/norm(pars.inversObsl.* pars.ox_signal);

iterate = 0;

end

