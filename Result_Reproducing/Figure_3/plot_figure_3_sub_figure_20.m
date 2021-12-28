function plot_figure_3_sub_figure_20(solver_choice, recovered_series)
% Function to generate the sub-figure in Figure 3, for each solver
% The coefficient reconstruction method and code script comes from paper:
%  L. Condat and A. Hirabayashi, Cadzow denoising upgraded: a new projection 
%    method for the recovery of Dirac pulse from noisy linear measurements, 
%    Samp. Theory Signal Image Processing, 14 (2015), 17–47.

if solver_choice == 1
    solver = 'FIHT';
elseif solver_choice == 2
    solver = 'PGD';
elseif solver_choice == 3
    solver = 'pMAP';
end

tau = 1;
M = 249;
N = length(recovered_series(:, 1));
r = 20;
figure_font_size = 20;

ox1 = recovered_series(:, 1);
recovered_signal = recovered_series(:, solver_choice + 1);

Hp0 = hankel(ox1(1:(N-r)), ox1((N-r):end));
[~, ~, V1] = svd(Hp0, 0);
tk=-angle(roots(flipud(V1(:,r+1))))*tau/2/pi;
tk = sort(tk);
ind=find(tk<0);
tk(ind)=tk(ind)+tau;
ak=real(pinv(exp(-1i*2*pi/tau*(-M:M)'*tk'))*ox1);

    f = figure('Name',sprintf('Solver %s, rank = 20', solver), 'NumberTitle','off', 'visible','off');

    vec=ak(1)*(1+2*sum(cos(2*pi*(1:M)'*((0:5000)/5000-tk(1)/tau)),1))/N;
    for k=2:length(tk)
        vec=vec+ak(k)*(1+2*sum(cos(2*pi*(1:M)'*((0:5000)/5000-tk(k)/tau)),1))/N;
    end
    plot((0:5000)/5000, vec,'color',[0.8,0.8,0.8]);
    hold on
    h1=stem(tk/tau,ak,'Color','k','MarkerEdgeColor','k','Markersize',7);
    xlim([0 1]);
    ymin=min(0,min(ak));
    ymax=max(0,max(ak));
    ylim([ymin-0.1 ymax+0.1]);
    set(get(h1,'BaseLine'),'LineStyle',':');
    set(gca,'FontSize',figure_font_size);
    xlabel('Frequency','fontsize',figure_font_size);
    ylabel('Amplitude','fontsize',figure_font_size)
    
    
Hp1 = hankel(recovered_signal(1:(N-r)), recovered_signal((N-r):end));

[~, ~, V2] = svd(Hp1, 0);
estimtk=-angle(roots(flipud(V2(:,r+1))))*tau/2/pi;
estimtk = sort(estimtk);

ind=find(estimtk<0);
estimtk(ind)=estimtk(ind)+tau;			
estimak=real(pinv(exp(-1i*2*pi/tau*(-M:M)'*estimtk'))*ox1);
edgecolor='r';

    ymin=min(ymin,min(estimak));
    ymax=max(ymax,max(estimak));
    h2 = stem(estimtk/tau,estimak,'--','Color',edgecolor,'MarkerFaceColor',edgecolor, ...
          'Marker'   , 'h' , 'MarkerEdgeColor',edgecolor,'Markersize',7);
    
SLegend = legend( ...
    [h1, h2], ...
    'True Location', solver, ...
    'location', 'northeast');
set(SLegend, 'FontSize', 16);

if solver_choice == 1
    out_put_file_name = 'Output/Fig_3(a).eps';
elseif solver_choice == 2
    out_put_file_name = 'Output/Fig_3(c).eps';
else
    out_put_file_name = 'Output/Fig_3(e).eps';
end

saveas(f,out_put_file_name,'epsc')
end

