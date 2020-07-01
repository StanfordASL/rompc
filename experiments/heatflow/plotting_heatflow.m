close all; clc; clear;
path = fileparts(which('plotting_heatflow.m'));
load(strcat(path, '/data/DATA_ROMPC.mat'));
figpath = strcat(path,'/figures/');
colors = {'b', [0.75, 0, 0.75], [0 0.6 0.3], 'r', 'k'};
markers = {'o','x','v','s','d'};
mark_idx = round(linspace(1,DATA_ROMPC.T-1,15));
r = DATA_ROMPC.opt.setpoint.r;

% Compute cost differences
fprintf('ROMPC simulation yielded a cost of %.2f.\n', DATA_ROMPC.J);


%% Performance Variable Response
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Heatflow','Interpreter','latex','Fontsize',22);
plot(DATA_ROMPC.z(1,:), 'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.z(2,:), 'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.z(3,:), 'color',colors{3},'marker',markers{3},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.z(4,:), 'color',colors{4},'marker',markers{4},'markerindices',mark_idx,'Linewidth',1);
plot(r(1)*ones(1,DATA_ROMPC.T), 'color',colors{1},'Linewidth',1,'Linestyle','--');
plot(r(2)*ones(1,DATA_ROMPC.T), 'color',colors{2},'Linewidth',1,'Linestyle','--');
plot(r(3)*ones(1,DATA_ROMPC.T), 'color',colors{3},'Linewidth',1,'Linestyle','--');
plot(r(4)*ones(1,DATA_ROMPC.T), 'color',colors{4},'Linewidth',1,'Linestyle','--');
ylim([-0.4, 1.2]);
xlabel('Time, $k$','Interpreter','latex','FontSize',22);
ylabel('Performance Output, $z(k)$','Interpreter','latex','FontSize',22);
legend({'$z_1$', '$z_2$', '$z_3$', '$z_4$'}, 'Interpreter','latex',...
                'FontSize',18,'Location','northeast','Orientation','horizontal');
legend('boxoff');            
% filename = strcat(figpath, 'z_heatflow');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')


%% Plot Control
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Heatflow','Interpreter','latex','Fontsize',22);
uUB = 100;
uLB = -100;
plot(DATA_ROMPC.u(1,:),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u(2,:),'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u(3,:),'color',colors{3},'marker',markers{3},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u(4,:),'color',colors{4},'marker',markers{4},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u(5,:),'color',colors{5},'marker',markers{5},'markerindices',mark_idx,'Linewidth',1);
plot(uUB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
plot(uLB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
ylim([-105, 130]);
xlabel('Time, $k$','Interpreter','latex','FontSize',22);
ylabel('Control Input, $u(k)$','Interpreter','latex','FontSize',22);
legend({'$u_1$', '$u_2$', '$u_3$', '$u_4$', '$u_5$', 'constraints'}, 'Interpreter','latex',...
            'FontSize',18,'Location','northeast','Orientation','horizontal');
legend('boxoff');
% filename = strcat(figpath, 'u_heatflow');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')



