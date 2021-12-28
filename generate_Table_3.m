%% Script to generate Table 3 in submitted paper
%  This table records the kpi performances of candidate solvers in taclking
%   spectral sparse signal when some observations are affected by noises.
%  Candidate solvers including Cadzow, DRI, FIHT, PGD and pMAP. 
%  For detailed information and references, please see README

clear all
clc
AddPath();
rng('default')
format short

%% Generate Table 3

instance_amount = 50;

kpi_table_3 = experiment_table_3(instance_amount);

%% Save result of table 3 into a csv file

fprintf('\n-------------------------------------------------');
fprintf('\n Recording numerical results in Table 3');

header_table = {'Cadzow' 'DRI' 'FIHT' 'PGD' 'pMAP'};
comma_header_table = [header_table;repmat({','},1,numel(header_table))];
comma_header_table = comma_header_table(:)';
text_header_table = cell2mat(comma_header_table); 

fid = fopen('Output/result_table_3.csv','w'); 
fprintf(fid,'%s\n',text_header_table);
fclose(fid);

dlmwrite('Output/result_table_3.csv', kpi_table_3, '-append');

fprintf('\n-------------------------------------------------');
fprintf('\n Table 3 recorded sucessfully, result can be seen from result_table_3.csv in /Output folder');
fprintf('\n-------------------------------------------------');