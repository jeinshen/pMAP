%% Script to generate Table 1 in submitted paper
%  This table records the kpi performances of candidate solvers in taclking
%   real valued time series denoising prolem
%  Candidate solvers including Cadzow, DRI and pMAP. 
%  For detailed information and references, please see README

clear all;
clc
AddPath();

rng('default')
format short

%% Generate results used in Table 1
instance_amount = 50;

kpi_table_1 = experiment_table_1(instance_amount);

%% Save result of table 1 into a csv file

fprintf('\n-------------------------------------------------');
fprintf('\n Recording numerical results in Table 1');

header_table_1 = {'pMAP' 'Cadzow' 'DRI'};
comma_header_table_1 = [header_table_1;repmat({','},1,numel(header_table_1))];
comma_header_table_1 = comma_header_table_1(:)';
text_header_table_1 = cell2mat(comma_header_table_1);

fid = fopen('Output/result_table_1.csv','w'); 
fprintf(fid,'%s\n',text_header_table_1);
fclose(fid);

dlmwrite('Output/result_table_1.csv', kpi_table_1, '-append');

fprintf('\n-------------------------------------------------');
fprintf('\n Table 1 recorded sucessfully, result can be seen from result_table_1.csv from /Output folder');
fprintf('\n-------------------------------------------------');
