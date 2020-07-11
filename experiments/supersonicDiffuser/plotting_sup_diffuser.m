close all; clc; clear;
path = fileparts(which('plotting_sup_diffuser.m'));
load(strcat(path, '/data/DATA_ROMPC.mat'));
load(strcat(path, '/data/DATA_LQR.mat'));
figpath = strcat(path,'/figures/');
colors = {'b', [0.75, 0, 0.75], [0 0.6 0.3], 'r', 'k'};
markers = {'o','x','v','s','d'};
mark_idx = round(linspace(1,DATA_ROMPC.T-1,15));
dt = 0.025;
t_ROMPC = linspace(0, dt*DATA_ROMPC.T, DATA_ROMPC.T);
t_LQR = linspace(0, dt*DATA_LQR.T, DATA_LQR.T);
rho0 = 1.225; % nominal density

% Compute cost differences
fprintf('ROMPC simulation yielded a cost of %.2f.\n', DATA_ROMPC.J);
fprintf('ROLQR simulation yielded a cost of %.2f.\n', DATA_LQR.J);


%% Performance Variable Response
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Supersonic Diffuser','Interpreter','latex','Fontsize',22);
plot(t_ROMPC, DATA_ROMPC.z(1,:), 'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(t_LQR, DATA_LQR.z(1,:), 'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
% ylim([]);
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Throat Mach Perturbation, $\Delta M_T$','Interpreter','latex','FontSize',22);
legend({'$\Delta M_T$ (ROMPC)','$\Delta M_T$ (ROLQR)'}, 'Interpreter','latex',...
                'FontSize',18,'Location','northeast','Orientation','vertical');
legend('boxoff');            
% filename = strcat(figpath, 'z_sup_diffuser');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')


%% Plot Control
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Supersonic Diffuser','Interpreter','latex','Fontsize',22);
uUB = 3;
uLB = -1;
plot(t_ROMPC, 100*DATA_ROMPC.u(1,:),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(t_LQR, 100*DATA_LQR.u(1,:),'color',colors{1},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
% plot(uUB*ones(1,DATA.T), 'color','k','Linewidth',1,'Linestyle','--');
plot(t_ROMPC, uLB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
ylim([-2, 1]);
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Inlet Bleed Variation, $\Delta b$ [\%$\dot{m}_0$]','Interpreter','latex','FontSize',22);
legend({'$\Delta b$ (ROMPC)', '$\Delta b$ (ROLQR)', 'constraint'}, 'Interpreter','latex',...
            'FontSize',18,'Location','northeast','Orientation','vertical');
legend('boxoff');
% filename = strcat(figpath, 'u_sup_diffuser');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')


%% Performance Variable Error
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Supersonic Diffuser','Interpreter','latex','Fontsize',22);
plot(t_ROMPC, DATA_ROMPC.z(1,:)-DATA_ROMPC.zbar(1,:), 'color',colors{2},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
% ylim([]);
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Tracking Error, $\delta (\Delta M_T)$','Interpreter','latex','FontSize',22);           
% filename = strcat(figpath, 'dz_sup_diffuser');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')

