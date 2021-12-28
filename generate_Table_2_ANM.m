%% Script to generate Table 2 in submitted paper
%  This table records the kpi performances of candidate solvers in taclking
%   noiseless spectral sparse signal.
%  This script is used for ANM solver only
%  For detailed information and references, please see README

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

kpi_table_2_1 = experiment_table_2_1_anm(instance_amount);

kpi_table_2_2 = experiment_table_2_2_anm(instance_amount);

%% Save result of table 2 (ANM only) into a csv file

fprintf('\n-------------------------------------------------');
fprintf('\n Recording numerical results in Table 2 (ANM only)');

header_table_1 = {'ANM'};
comma_header_table_1 = [header_table_1;repmat({','},1,numel(header_table_1))];
comma_header_table_1 = comma_header_table_1(:)';
text_header_table_1 = cell2mat(comma_header_table_1);

fid = fopen('Output/result_table_2_1_anm.csv','w'); 
fprintf(fid,'%s\n',text_header_table_1);
fclose(fid);

dlmwrite('Output/result_table_2_1_anm.csv', kpi_table_2_1, '-append');

fid = fopen('Output/result_table_2_2_anm.csv','w'); 
fprintf(fid,'%s\n',text_header_table_1);
fclose(fid);

dlmwrite('Output/result_table_2_2_anm.csv', kpi_table_2_2, '-append');

fprintf('\n-------------------------------------------------');
fprintf('\n Table 2 (ANM only) recorded sucessfully, result can be seen from result_table_2_1_anm.csv and result_table_2_2_anm.csv from /Output folder');
fprintf('\n-------------------------------------------------');
