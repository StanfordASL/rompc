close all; clc; clear;
path = fileparts(which('plotting_small_synthetic.m'));
load(strcat(path, '/data/DATA_ROMPC.mat'));
load(strcat(path, '/data/DATA_LQR.mat'));
figpath = strcat(path,'/figures/');
colors = {'b', [0.75, 0, 0.75], [0 0.6 0.3], 'r', 'k'};
markers = {'o','x','v','s','d'};
mark_idx = round(linspace(1,DATA_ROMPC.T-1,15));

% Compute cost differences
fprintf('ROMPC simulation yielded a cost of %.2f.\n', DATA_ROMPC.J);


%% Performance Variable Response
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Small Synthetic','Interpreter','latex','Fontsize',22);
zUB = -50;
plot(DATA_ROMPC.z(1,:), 'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.z(2,:), 'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_LQR.z(1,:), 'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
plot(DATA_LQR.z(2,:), 'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
plot(zUB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
ylim([-90, 20]);
xlabel('Time, $k$','Interpreter','latex','FontSize',22);
ylabel('Performance Output, $z(k)$','Interpreter','latex','FontSize',22);
legend({'$z_1$ (ROMPC)', '$z_2$ (ROMPC)','$z_1$ (ROLQR)', '$z_2$ (ROLQR)','constraint'}, 'Interpreter','latex',...
                'FontSize',18,'Location','southeast','Orientation','vertical');
legend('boxoff');            
% filename = strcat(figpath, 'z_small_synthetic');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')


%% Plot Control
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Small Synthetic','Interpreter','latex','Fontsize',22);
uUB = 20;
uLB = -20;
plot(DATA_ROMPC.u(1,:),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_LQR.u(1,:),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
% plot(uUB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
plot(uLB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
% ylim([]);
xlabel('Time, $k$','Interpreter','latex','FontSize',22);
ylabel('Control Input, $u(k)$','Interpreter','latex','FontSize',22);
legend({'$u_1$ (ROMPC)', '$u_1$ (ROLQR)', 'constraints'}, 'Interpreter','latex',...
            'FontSize',18,'Location','southeast','Orientation','vertical');
legend('boxoff');
% filename = strcat(figpath, 'u_small_synthetic');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')



