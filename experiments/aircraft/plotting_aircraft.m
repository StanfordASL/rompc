close all; clc; clear;
path = fileparts(which('plotting_aircraft.m'));
load(strcat(path, '/data/DATA.mat'));
figpath = strcat(path,'/figures/');
colors = {'b', [0.75, 0, 0.75], [0 0.6 0.3], 'r', 'k'};
markers = {'o','x','v','s','d'};
T = 3000;
mark_idx = round(linspace(1,T-1,15));


%% Performance Variable Response
figure('color',[1,1,1],'Position', [1, 1, 2000,500]); hold on;
subplot(1,4,1); hold on;
title('Position','Interpreter','latex','Fontsize',22);
plot(DATA.t(1:T), DATA.z(1,1:T), 'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA.t(1:T), DATA.z(2,1:T), 'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
% ylim([]);
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Position Error [m]','Interpreter','latex','FontSize',22);
legend({'$\delta p_x$', '$\delta p_z$'}, 'Interpreter','latex',...
                'FontSize',18,'Location','northeast','Orientation','vertical');
legend('boxoff');

subplot(1,4,2); hold on;
title('Velocity','Interpreter','latex','Fontsize',22);
plot(DATA.t(1:T), DATA.z(4,1:T),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA.t(1:T), DATA.z(5,1:T),'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Velocity Error, (m)','Interpreter','latex','FontSize',22);
legend({'$\delta \dot{p}_x$', '$\delta \dot{p}_z$'}, 'Interpreter','latex',...
                'FontSize',18,'Location','best','Orientation','vertical');
legend('boxoff');       

subplot(1,4,3); hold on;
title('Pitch','Interpreter','latex','Fontsize',22);
plot(DATA.t(1:T), rad2deg(DATA.z(3,1:T)),'color',colors{3},'marker',markers{3},'markerindices',mark_idx,'Linewidth',1);
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Pitch Error, $(^\circ)$','Interpreter','latex','FontSize',22);
legend({'$\delta \theta$'}, 'Interpreter','latex',...
                'FontSize',18,'Location','best','Orientation','vertical');
legend('boxoff');       

subplot(1,4,4); hold on;
title('Pitch Rate','Interpreter','latex','Fontsize',22);
plot(DATA.t(1:T), rad2deg(DATA.z(6,1:T)),'color',colors{3},'marker',markers{4},'markerindices',mark_idx,'Linewidth',1);
plot(DATA.t(1:T), 10*ones(1,T), 'color','k','Linewidth',1,'Linestyle','--');
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Pitch Rate Error, ($^\circ$/s)','Interpreter','latex','FontSize',22);
ylim([-2, 11]);
legend({'$\delta \dot{\theta}$','constraint'}, 'Interpreter','latex',...
                'FontSize',18,'Location','east','Orientation','vertical');
legend('boxoff');
% filename = strcat(figpath, 'z_aircraft');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')


%% Plot Control
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Aircraft','Interpreter','latex','Fontsize',22);
plot(DATA.t(1:T),DATA.u(1,1:T),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA.t(1:T),DATA.u(2,1:T),'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
% ylim([]);
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Control Inputs, $u(t)$','Interpreter','latex','FontSize',22);
legend({'$T$, [N]', '$M$, [Nm]'}, 'Interpreter','latex',...
            'FontSize',18,'Location','northeast','Orientation','vertical');
legend('boxoff');
% filename = strcat(figpath, 'u_aircraft');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')



