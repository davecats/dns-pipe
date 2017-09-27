%% Clear Command Window and Workspace
%  -------------------------------------------------------------------
clc
clear variables
restoredefaultpath

%% Define what to plot
printplots=0;
%% Add subforlder ./../base/ into MATLAB search paths
%  -------------------------------------------------------------------
addpath('./../base/');

%% Give the path 
% filepath='../../DATA/nx-32-ny-128-nz-128-alfa0-alpha0-2.0/';
% filepath='../../DATA/Re5000-nx-128-ny-256-nz-256-alpha0-2.0/';
filepath_short='Re28000/nx-256-ny-500-nz-512-alfa0-4.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-0.5/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-2.0/';
% filepath_short='Re5000-nx-128-ny-256-nz-256-alpha0-2.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-2/';
% filepath_short='Re5000-nx-128-ny-256-nz-256-alpha0-2.0/';
load(['../../../Results/',filepath_short,'Statistics.mat'],'stats');


%% Define Plots
sizeoffonts=45;
widthlines=3;
colorlist={[0,150,130]/255,[162 34 35]/255,[70 100 170]/256,[252 229 0]/255,[140 182 60]/256,[223 155 27]/255,[167 130 46]/255,[163 16 124]/255,[35 161 224]/255};
set(groot, 'defaultAxesTickLabelInterpreter','latex')



    fig3=figure;
    set(fig3,'Position',[10 50 1050 800],'Color','w')
    set(fig3, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
    ax3 = axes('Parent', fig3);
%     pax3_1=plot(ax3,stats.y,stats.y,'o-','LineWidth',widthlines,'Color',colorlist{1});
    xlim([0,1])
    ylim([0,1.1])
    hold on
    set(ax3,'FontSize',sizeoffonts)
    xlabel('$(r-r_i)/\delta_{99}$ [-]','interpreter','latex','fontsize',sizeoffonts);
    ylabel('$u_{\varphi}/u_{\varphi}(r_i)$ [-]','interpreter','latex','fontsize',sizeoffonts);
    box on


%%
counter=0;
colorcount=0;
for iF=1:numel(stats.t)%numel(filenamelist)-1
    iF
    %Check colorcount
    if colorcount<9;colorcount=colorcount+1; else; colorcount=1; end

        plot(ax3,(stats.y(stats.d99_pos(iF):numel(stats.y))-stats.d99(iF))/(1-stats.d99(iF)),1/stats.d99(iF)*squeeze(real(stats.Vm(3,stats.d99_pos(iF):end,iF))),'-','LineWidth',widthlines,'Color',colorlist{colorcount});
         
        if printplots
            mkdir(['../../../Results/',filepath_short,'normal_ruben2']);
            print(fig3,'-dpng',['../../../Results/',filepath_short,'normal_ruben2/',num2str(10000+1000*stats.t(iF)),'.png'],'-r100');
        end
  
end
