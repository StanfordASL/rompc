close all; clc; clear;
path = fileparts(which('plotting_dist_column.m'));
load(strcat(path, '/data/DATA_ROMPC.mat'));
load(strcat(path, '/data/DATA_LQR.mat'));
figpath = strcat(path,'/figures/');
colors = {'b', [0.75, 0, 0.75], [0 0.6 0.3], 'r', 'k'};
markers = {'o','x','v','s','d'};
mark_idx = round(linspace(1,DATA_ROMPC.T-1,15));
dt = 1;

% Compute cost differences
fprintf('ROMPC simulation yielded a cost of %.2f.\n', DATA_ROMPC.J);
fprintf('ROLQR simulation yielded a cost of %.2f.\n', DATA_ROMPC.J);

%% Performance Variable Response
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Distillation Column','Interpreter','latex','Fontsize',22);
plot(DATA_ROMPC.z(1,:), 'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.z(2,:), 'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_LQR.z(1,:), 'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
plot(DATA_LQR.z(2,:), 'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
plot(-.01*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
% ylim([]);
xlabel('Time, [min]','Interpreter','latex','FontSize',22);
ylabel('Product Compositions, $y_D$, $x_B$ [mol/mol]','Interpreter','latex','FontSize',22);
legend({'$y_D$ (ROMPC)', '$x_B$ (ROMPC)', '$y_D$ (ROLQR)', '$x_B$ (ROLQR)','constraint'}, 'Interpreter','latex',...
                'FontSize',18,'Location','east','Orientation','vertical');
legend('boxoff');            
% filename = strcat(figpath, 'z_dist_column');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')

%% Plot Control
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Distillation Column','Interpreter','latex','Fontsize',22);
uUB = 1;
uLB = -1;
plot(DATA_ROMPC.z(5,:),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.z(6,:),'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.z(7,:),'color',colors{3},'marker',markers{3},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.z(8,:),'color',colors{4},'marker',markers{4},'markerindices',mark_idx,'Linewidth',1);
plot(uUB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
plot(uLB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
ylim([-1.05, 1.05]);
xlim([0, 50]);
xlabel('Time, [min]','Interpreter','latex','FontSize',22);
ylabel('Control Inputs, $u$ [kmol/min]','Interpreter','latex','FontSize',22);
legend({'$L$', '$V$', '$D$', '$B$', 'constraints'}, 'Interpreter','latex',...
            'FontSize',18,'Location','northeast','Orientation','horizontal');
legend('boxoff');
% filename = strcat(figpath, 'u_dist_column');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')


%% Plot Control
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Distillation Column','Interpreter','latex','Fontsize',22);
uUB = 0.2;
uLB = -0.2;
plot(DATA_ROMPC.u(1,:),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u(2,:),'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u(3,:),'color',colors{3},'marker',markers{3},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u(4,:),'color',colors{4},'marker',markers{4},'markerindices',mark_idx,'Linewidth',1);
plot(uUB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
plot(uLB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
ylim([-0.21, 0.25]);
xlim([0, 50]);
xlabel('Time, [min]','Interpreter','latex','FontSize',22);
ylabel('Control Rates, $\Delta u$ [kmol/min$^2$]','Interpreter','latex','FontSize',22);
legend({'$\Delta L$', '$\Delta V$', '$\Delta D$', '$\Delta B$', 'constraints'}, 'Interpreter','latex',...
            'FontSize',18,'Location','northeast','Orientation','horizontal');
legend('boxoff');
% filename = strcat(figpath, 'du_dist_column');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')



