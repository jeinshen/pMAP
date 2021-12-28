%% Script to generate Table 2 in submitted paper
%  This table records the kpi performances of candidate solvers in taclking
%   noiseless spectral sparse signal.

clear all
clc
AddPath();
rng('default')
format short

%% Generate Table 2
%   kpi_table_2_1 contains the first half of table 2 where 30% data is
%       known
%   kpi_table_2_2 contains the second half of table 2 where 60% data is
%       known

instance_amount = 50;

kpi_table_2_1 = experiment_table_2_1(instance_amount);

kpi_table_2_2 = experiment_table_2_2(instance_amount);

%% Save result of table 2 (without ANM) into a csv file

fprintf('\n-------------------------------------------------');
fprintf('\n Recording numerical results in Table 2 (without ANM)');

header_table = {'Cadzow' 'DRI' 'FIHT' 'PGD' 'pMAP'}; 
comma_header_table = [header_table;repmat({','},1,numel(header_table))]; 
comma_header_table = comma_header_table(:)';
text_header_table = cell2mat(comma_header_table); 

fid = fopen('Output/result_table_2_1.csv','w'); 
fprintf(fid,'%s\n',text_header_table);
fclose(fid);

dlmwrite('Output/result_table_2_1.csv', kpi_table_2_1, '-append');

fid = fopen('Output/result_table_2_2.csv','w'); 
fprintf(fid,'%s\n',text_header_table);
fclose(fid);

dlmwrite('Output/result_table_2_2.csv', kpi_table_2_2, '-append');

fprintf('\n-------------------------------------------------');
fprintf('\n Table 2 (without ANM) recorded sucessfully, result can be seen from result_table_2_1.csv and result_table_2_2.csv in /Output folder');
fprintf('\n-------------------------------------------------');
