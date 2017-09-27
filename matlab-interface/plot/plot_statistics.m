%% This function plots different statistics calculated by statistics/postpro_statistics.m


%% Clear Command Window and Workspace
%  -------------------------------------------------------------------
clc
clear variables
restoredefaultpath
%% Define what to plot

%plots the Reynoldsstress uphiuphi, uphiur in wall units
plotvar_t=0;
%loglayer: plots the velocity profile in wallunits
plotlog=0;
% normalizes the velocity profile like in the paper of verschoof et al.
ruben=1;
%plots the Reynoldsstress uphiuphi in wall units
varphi=0;
% plots u_phi and the Reynoldsstresses and the boundary layer thickness
profiles=0;
% plots omega and uphi over r
plotomega=0;

%% If this is activated the plots will be saved as pngs
printplots=1;
% Maximum position of the reynolds stress uphiuphi (used in some plots) 
maxyplus=6;
%% Add subforlder ./base/ into MATLAB search paths
%  -------------------------------------------------------------------
addpath('./../base/');

% filepath='../../DATA/nx-32-ny-128-nz-128-alfa0-alpha0-2.0/';
% filepath='../../DATA/Re5000-nx-128-ny-256-nz-256-alpha0-2.0/';
filepath_short='Re28000/nx-256-ny-500-nz-512-alfa0-4.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-0.5/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-2.0/';
% filepath_short='Re5000-nx-128-ny-256-nz-256-alpha0-2.0/';
% filepath_short='Re3000/nx-32-ny-128-nz-128-alfa0-alpha0-8.0-omega-2/';
% filepath_short='Re5000-nx-128-ny-256-nz-256-alpha0-2.0/';
load(['./../../../Results/',filepath_short,'Statistics.mat'],'stats');


%% Define Plots
sizeoffonts=30;
widthlines=3;
colorlist={[0,150,130]/255,[162 34 35]/255,[70 100 170]/256,[252 229 0]/255,[140 182 60]/256,[223 155 27]/255,[167 130 46]/255,[163 16 124]/255,[35 161 224]/255};
set(groot, 'defaultAxesTickLabelInterpreter','latex')

if plotvar_t
    % Variance uphiuphi, uphiur
    fig1=figure;
    set(fig1,'Position',[10 50 800 800],'Color','w')
    set(fig1, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
    ax1 = axes('Parent', fig1);
    pax1_1=plot(ax1,stats.y,stats.y,'o-','LineWidth',widthlines,'Color',colorlist{1});
    hold on
    pax1_2=plot(ax1,(maxyplus)*ones(11,1),0:1:10,'LineWidth',widthlines,'Color',colorlist{2});
    pax1_3=plot(ax1,1-stats.y,stats.y,'o-','LineWidth',widthlines,'Color',colorlist{3});
    xlim([10^-1,10^3])
    maxfig1=max(max(squeeze(stats.var(3,3,:,:))/abs(stats.tau_w(:)')));
    ylim([-0.1,1.3*maxfig1])
    xticks([1 10 100])
    grid on
    set(ax1,'xscale','log')
    set(ax1,'FontSize',sizeoffonts)
    legend([pax1_1,pax1_3,pax1_2],{'$<u_{\varphi}''u_{\varphi}''>/u_{\tau}^2$',...
        '$<u_{\varphi}''u_{r}''>/u_{\tau}^2$','max$(<u_{\varphi}''u_{\varphi}''>)$'}, ...
        'interpreter','latex','fontsize',sizeoffonts,'Location','northeast')
    xlabel('$y+$ [-]','interpreter','latex','fontsize',sizeoffonts);
end


% Loglayer
if plotlog
    fig2=figure;
    set(fig2,'Position',[10 50 800 800],'Color','w')
    set(fig2, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
    ax2 = axes('Parent', fig2);
    pax2_1=plot(ax2,stats.y,stats.y,'o-','LineWidth',widthlines,'Color',colorlist{1});
    hold on
    pax2_2=plot(ax2,stats.y,stats.y,'LineWidth',widthlines,'Color',colorlist{2});
    ylim([0,20])
    xlim([10^-1,10^3])
    xticks([1 10 100])
    set(ax2,'xscale','log')
    set(ax2,'FontSize',sizeoffonts)
    legend([pax2_1,pax2_2],{'$u_{\varphi}/u_{\tau}$','subl.+log'},'interpreter','latex','fontsize',sizeoffonts,'Location','northwest')
    xlabel('$y+$ [-]','interpreter','latex','fontsize',sizeoffonts);
    % ylabel('$u_{\varphi}/u_{\tau}$ [-]','interpreter','latex','fontsize',sizeoffonts);
end
% Normalized velocity fields
if ruben
    fig3=figure;
    set(fig3,'Position',[10 50 800 800],'Color','w')
    set(fig3, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
    ax3 = axes('Parent', fig3);
    pax3_1=plot(ax3,stats.y,stats.y,'o-','LineWidth',widthlines,'Color',colorlist{1});
    xlim([0,1])
    ylim([0,1.1])
    set(ax3,'FontSize',sizeoffonts)
    xlabel('$(r-r_i)/(R-r_i)$ [-]','interpreter','latex','fontsize',sizeoffonts);
    ylabel('$u_{\varphi}(r)/u_{\varphi}(r_i)$ [-]','interpreter','latex','fontsize',sizeoffonts);
end

if varphi
    % Variance uphiuphi one plot
    fig4=figure;
    set(fig4,'Position',[10 50 800 800],'Color','w')
    set(fig4, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
    ax4 = axes('Parent', fig4);
    hold on
    xlim([10^-1,10^3])
    maxfig4=max(max(squeeze(stats.var(3,3,:,:))'./abs(stats.tau_w(:))));
    ylim([-0.1,1.3*maxfig4])
    set(ax4,'xscale','log')
    grid on
    xticks([1 10 100])
    set(ax4,'FontSize',sizeoffonts)
    xlabel('$y+$ [-]','interpreter','latex','fontsize',sizeoffonts);
    ylabel('$<u_{\varphi}''u_{\varphi}''>/u_{\tau}^2$','interpreter','latex','fontsize',sizeoffonts);
    pax4_2=plot(ax4,(maxyplus)*ones(11,1),0:1:10,'LineWidth',widthlines,'Color',colorlist{2});
end
% Variance uphiuphi one plot
if profiles
    fig5=figure;
    set(fig5,'Position',[10 50 800 800],'Color','w')
    set(fig5, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
    ax5 = axes('Parent', fig5);
    hold on
    xlim([0,1])
    ylim([-0.1,1.2])
    set(ax5,'FontSize',sizeoffonts)
    xlabel('$r/R$ [-]','interpreter','latex','fontsize',sizeoffonts);
    pax5_1=plot(ax5,stats.y,squeeze(stats.Vm(3,:,1)),'LineWidth',widthlines,'Color',colorlist{1});
    hold on
    grid on
    norm1=max(max(squeeze(stats.var(3,3,:,:))));
    pax5_bl=plot(ax5,(1)*ones(11,1),0:0.1:1,'LineWidth',widthlines,'Color',colorlist{2});
    pax5_2=plot(ax5,stats.y,stats.y,'LineWidth',widthlines,'Color',colorlist{3});
    pax5_3=plot(ax5,stats.y,stats.y,'LineWidth',widthlines,'Color',colorlist{4});
    pax5_4=plot(ax5,stats.y,stats.y,'LineWidth',widthlines,'Color',colorlist{5});
    pax5_5=plot(ax5,stats.y,stats.y,'LineWidth',widthlines,'Color',colorlist{6});
    pax5_6=plot(ax5,stats.y,stats.y,'Color','k');
    legend([pax5_1,pax5_2,pax5_3,pax5_4,pax5_5,pax5_bl],{'$<u_{\varphi}>\Omega R$', ...
        '$<u_{\varphi}''u_{\varphi}''>$', ...
        '$<u_{\varphi}''u_{r}''>$', ...
        '$<u_{z}''u_{z}''>$', ...
        '$<u_{r}''u_{r}''>$', ...
        '$\delta_{99}/R$'}, ...
        'interpreter','latex','fontsize',sizeoffonts,'Location','northwest')
end

if plotomega
    fig6=figure;
    set(fig6,'Position',[10 50 800 800],'Color','w')
    set(fig6, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
    ax6 = axes('Parent', fig6);
    hold on
    xlim([0,1])
    ylim([-0.5,1.1])
    set(ax6,'FontSize',sizeoffonts)
    xlabel('$r/R$ [-]','interpreter','latex','fontsize',sizeoffonts);
    pax6_1=plot(ax6,stats.y,squeeze(stats.Vm(3,:,1)),'LineWidth',widthlines,'Color',colorlist{1});
    hold on
    grid on
    norm1=max(max(squeeze(stats.var(3,3,:,:))));
    pax6_bl=plot(ax6,(1)*ones(11,1),-1:0.2:1,'LineWidth',widthlines,'Color',colorlist{2});
    pax6_2=plot(ax6,stats.y(3:end),stats.y(3:end),'LineWidth',widthlines,'Color',colorlist{3});
    
    legend([pax6_1,pax6_2,pax6_bl],{'$<u_{\varphi}>/\Omega R$', ...
        '$<\omega_z>/2 \Omega$', ...
        '$\delta_{99}/R$'}, ...
        'interpreter','latex','fontsize',sizeoffonts,'Location','southwest')
end

%%
% fig=figure;
counter=0;
colorcount=0;
for iF=1500:numel(stats.t)%numel(filenamelist)-1
    iF
    %Check colorcount
    if colorcount<9;colorcount=colorcount+1; else; colorcount=1; end
    %fig1
    if plotvar_t
        set(pax1_1,'Xdata',stats.yplus(:,iF),'Ydata',squeeze(stats.var(3,3,:,iF))/abs(stats.tau_w(iF)));
        set(pax1_3,'Xdata',stats.yplus(:,iF),'Ydata',squeeze(stats.var(2,3,:,iF))/abs(stats.tau_w(iF)));
        mkdir(['../../Results/',filepath_short,'var_phiphi_phir']);
        print(fig1,'-dpng',['../../Results/',filepath_short,'var_phiphi_phir/',num2str(10000+1000*stats.t(iF)),'.png'],'-r100');
    end
    if plotlog
        %fig2: Loglayer
        numelsublayer=numel(find(stats.yplus(:,iF)<10));
        uplus=zeros(size(stats.yplus(:,iF)));
        uplus(end-numelsublayer+1:end)=stats.yplus(end-numelsublayer+1:end,iF);
        kappa=0.41; cplus=5;
        uplus(1:end-numelsublayer)=1/kappa*log(stats.yplus(1:end-numelsublayer,iF))+cplus;
        set(pax2_1,'Xdata',stats.yplus(:,iF),'Ydata',squeeze(real(stats.Vm(3,:,iF)))/abs(stats.u_tau(iF)));
        set(pax2_2,'Xdata',stats.yplus(:,iF),'Ydata',uplus);
        if printplots
            mkdir(['../../../Results/',filepath_short,'loglayer']);
            print(fig2,'-dpng',['../../../Results/',filepath_short,'loglayer/',num2str(10000+1000*stats.t(iF)),'.png'],'-r100');
        end
    end
    if ruben
        %fig 3
        set(pax3_1,'Xdata',(stats.y(stats.d99_pos(iF):numel(stats.y))-stats.d99(iF))/(1-stats.d99(iF)),'Ydata',1/stats.d99(iF)*squeeze(real(stats.Vm(3,stats.d99_pos(iF):end,iF))));
        if printplots
            mkdir(['../../../Results/',filepath_short,'normal_ruben']);
            print(fig3,'-dpng',['../../../Results/',filepath_short,'normal_ruben/',num2str(10000+1000*stats.t(iF)),'.png'],'-r100');
        end
    end
    if varphi
        %fig 4
        plot(ax4,stats.yplus(:,iF),squeeze(stats.var(3,3,:,iF))/abs(stats.tau_w(iF)),'o-','LineWidth',widthlines,'Color',colorlist{colorcount});
    end
    if profiles
        %fig 5
        norm2=norm1; %max(squeeze(field.var(3,3,1,:)));
        set(pax5_1,'Ydata',real(stats.Vm(3,:,iF)));
        set(pax5_bl,'Xdata',(stats.d99(iF))*ones(11,1));
        set(pax5_2,'Ydata',squeeze(stats.var(3,3,:,iF))/norm2);
        set(pax5_3,'Ydata',squeeze(stats.var(2,3,:,iF))/norm2);
        set(pax5_4,'Ydata',squeeze(stats.var(1,1,:,iF))/norm2);
        set(pax5_5,'Ydata',squeeze(stats.var(2,2,:,iF))/norm2);
        if printplots
            mkdir(['../../../Results/',filepath_short,'/profiles']);
            print(fig5,'-dpng',['../../../Results/',filepath_short,'profiles/',num2str(10000+1000*stats.t(iF)),'.png'],'-r100');
        end
    end
    
    
    if plotomega
        set(pax6_1,'Ydata',real(stats.Vm(3,:,iF)));
        set(pax6_bl,'Xdata',(stats.d99(iF))*ones(11,1));
        set(pax6_2,'Ydata',real(stats.omegam(1,3:end,iF)/2));
        if printplots
            mkdir(['../../../Results/',filepath_short,'/omega']);
            print(fig6,'-dpng',['../../../Results/',filepath_short,'omega/',num2str(10000+1000*stats.t(iF)),'.png'],'-r100');
        end
    end
    
end
if varphi
    mkdir(['../../../Results/',filepath_short,'var_uphiphi']);
    print(fig4,'-dpng',['../../../Results/',filepath_short,'var_uphiphi/',num2str(10000+stats.t(iF-1)),'.png'],'-r100');
end
fig6=figure;
set(fig6,'Position',[10 50 800 800],'Color','w')
set(fig6, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
ax6 = axes('Parent', fig6);
pax6_1=plot(ax6,stats.t,stats.Re_tau,'o-','LineWidth',widthlines,'Color',colorlist{1});
hold on
set(ax6,'FontSize',sizeoffonts)
pax6_2=plot(ax6,stats.t,stats.Re_TKE,'o-','LineWidth',widthlines,'Color',colorlist{2});
pax6_3=plot(ax6,stats.t,stats.Re_TC,'o-','LineWidth',widthlines,'Color',colorlist{3});
xlabel('$\Omega t$ [-]','interpreter','latex','fontsize',sizeoffonts);
ylabel('$Re$ [-]','interpreter','latex','fontsize',sizeoffonts);
legend([pax6_1,pax6_2,pax6_3],{'$Re_{\tau}$','$Re_{TKE}$','$Re_{TC}$'}, ...
    'interpreter','latex','fontsize',sizeoffonts,'Location','northwest')
if printplots
    print(fig6,'-dpng',['../../../Results/',filepath_short,'Reynoldsnumbers.png'],'-r100');
    savefig(fig6,['../../../Results/',filepath_short,'Reynoldsnumbers.fig'])
end
fig7=figure;
set(fig7,'Position',[10 50 800 800],'Color','w')
set(fig7, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
ax7= axes('Parent', fig7);
pax7_1=plot(stats.t,squeeze(stats.ypluspos(3,3,:)),'o-','LineWidth',widthlines,'Color',colorlist{1});
xlabel('$\Omega t$ [-]','interpreter','latex','fontsize',sizeoffonts);
ylabel('$y^+$ [-]','interpreter','latex','fontsize',sizeoffonts);
legend(pax7_1,{'$y^+(max(<u_{\varphi}''u_{\varphi}''>))$'}, ...
    'interpreter','latex','fontsize',sizeoffonts,'Location','northwest');
set(ax7,'FontSize',sizeoffonts)
ylim([0.5*min(stats.ypluspos(3,3,11:end)),1.4*max(stats.ypluspos(3,3,11:end))])
if printplots
    print(fig7,'-dpng',['../../../Results/',filepath_short,'Varphiphi.png'],'-r100');
    savefig(fig7,['../../../Results/',filepath_short,'Varphiphi.fig'])
end
fig8=figure;
set(fig8,'Position',[10 50 800 800],'Color','w')
set(fig8, 'PaperUnits', 'inches', 'PaperPosition', [0.5, 0.2, 18, 3.00],'PaperPositionMode','auto');
ax8= axes('Parent', fig8);
pax8_1=plot(ax8,stats.t,squeeze(stats.ypluspos(3,3,:)),'o-','LineWidth',widthlines,'Color',colorlist{1});
hold on
set(ax8,'FontSize',sizeoffonts)
pax8_2=plot(ax8,stats.t,squeeze(stats.ypluspos(2,2,:)),'o-','LineWidth',widthlines,'Color',colorlist{2});
pax8_3=plot(ax8,stats.t,squeeze(stats.ypluspos(1,1,:)),'o-','LineWidth',widthlines,'Color',colorlist{3});
pax8_4=plot(ax8,stats.t,squeeze(stats.ypluspos(2,3,:)),'o-','LineWidth',widthlines,'Color',colorlist{4});
xlabel('$\Omega t$ [-]','interpreter','latex','fontsize',sizeoffonts);
ylabel('$y^+$ [-]','interpreter','latex','fontsize',sizeoffonts);
legend([pax8_1,pax8_2,pax8_3,pax8_4],{'$y^+(max(<u_{\varphi}''u_{\varphi}''>))$', ...
    '$y^+(max(<u_{r}''u_{r}''>))$', ...
    '$y^+(max(<u_{z}''u_{z}''>))$',...
    '$y^+(max(<u_{\varphi}''u_{r}''>))$'}, ...
    'interpreter','latex','fontsize',sizeoffonts,'Location','northeast');
xlim([0.3,stats.t(end)])
ylim([0.4*min(stats.ypluspos(3,3,:)),2*max(squeeze(stats.ypluspos(2,2,35:end))) ])
grid on
if printplots
    print(fig8,'-dpng',['../../../Results/',filepath_short,'Var_t3.png'],'-r100');
    savefig(fig8,['../../../Results/',filepath_short,'Var_t3.fig'])
    
end

