function figure_generator_figure_5_a(xerror_ratios)
% Function to plot figure 5.a
% This figure plot the relative gap between x^v and x^{v+1} at each iterate
%   for incomplete signal recovery experiment when input rank is
%   mis-approximated

fprintf('\n-------------------------------------------------');
fprintf('\n Step 3/4: generate Fig_5(a), approximation error plotting');


ratio_FIHT = xerror_ratios{1};
ratio_PGD = xerror_ratios{2};
ratio_pMAP = xerror_ratios{3};

f = figure('visible','off');
log_ratio_FIHT = log(ratio_FIHT);
iterations_fiht = max(size(ratio_FIHT));
iterations = 50;
font_size =24;

Line_FIHT = line((1 : iterations_fiht)', log_ratio_FIHT(1:iterations_fiht, 1)); 
hold on

log_ratio_PGD = log(ratio_PGD);
Line_PGD = line((1 : iterations)', log_ratio_PGD(1:iterations, 1)); 

log_ratio_pMAP = log(ratio_pMAP);
Line_pMAP = line((1 : iterations)', log_ratio_pMAP(1:iterations, 1)); 

set(Line_FIHT, ...                        
    'color'             ,  'r'    , ...
    'LineStyle'         ,  '--'         , ...
    'LineWidth'         , 1     , ...
    'Marker'            ,  's'          , ...
    'MarkerSize'        ,  8            , ...
    'MarkerEdgeColor'   ,  [.6 .6 .6]   , ...
    'MarkerFaceColor'   ,  [.7 .7 .7]);

set(Line_PGD, ...                         
    'color'             ,  'b'    , ...
    'LineStyle'         ,  ':'         , ...
    'LineWidth'         , 1     , ...
    'Marker'            ,  '.'          , ...
    'MarkerSize'        ,  20            , ...
    'MarkerEdgeColor'   ,  [.1 .1 1]   , ...
    'MarkerFaceColor'   ,  [.2 .3 .8]);
set(Line_pMAP, ...
    'color'             ,  [0.5 0 0.5]    , ...
    'LineStyle'         ,  '--'         , ...
    'LineWidth'         , 1     , ...
    'Marker'            ,  'd'          , ...
    'MarkerSize'        ,  8            , ...
    'MarkerEdgeColor'   ,  [.7 .1 .7]   , ...
    'MarkerFaceColor'   ,  [1 .1 .1]);
set(gca,'FontSize',font_size);
set(gca,'XTick',0:5:200);
ylabel('$log(\frac{\|x^{\nu+1} - x^{\nu}\|}{\|x^\nu\|})$', 'Interpreter','latex', 'fontsize',font_size);
xlabel('Iteration Step', 'fontsize',font_size);
set(gcf,'Position',[100 100 800 600]);
title('$log(\frac{\|x^{\nu+1} - x^{\nu}\|}{\|x^\nu\|})$ at each iteration', 'fontsize',font_size,  'Interpreter','latex')
SLegend = legend( ...
    [Line_FIHT, Line_PGD, Line_pMAP], ...
    'FIHT', 'PGD', 'pMAP',...
    'location', 'Best');
set(SLegend,'Interpreter','latex','FontSize', font_size);

saveas(f,'Output/Fig_5(a).eps','epsc')

fprintf('\n-------------------------------------------------');
fprintf('\n Figure 5.a succesfully generated. Please find the Fig_5(a).eps in /Output folder');
fprintf('\n-------------------------------------------------\n');


end

