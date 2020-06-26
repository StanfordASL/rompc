close all; clc; clear;
path = fileparts(which('plotting_tubular_reactor.m'));
load(strcat(path, '/data/DATA_ROMPC.mat'));
load(strcat(path, '/data/DATA_LQR.mat'));
figpath = strcat(path,'/figures/');
colors = {'b', [0.75, 0, 0.75], [0 0.6 0.3], 'r', 'k'};
markers = {'o','x','v','s','d'};
mark_idx = round(linspace(1,DATA_ROMPC.T-1,15));

% Compute cost differences
fprintf('ROMPC simulation yielded a cost of %.2f.\n', DATA_ROMPC.J);
fprintf('ROLQR simulation yielded a cost of %.2f.\n', DATA_LQR.J);

% Dimensionalize for better plotting
[FOM, ~, ~, ~, ~, ~, ~, PARAMS] = tubularReactor(false);
nf = size(FOM.Af,1);
DATA_ROMPC.u_dim = zeros(size(DATA_ROMPC.u));
for i = 1:DATA_ROMPC.T
    DATA_ROMPC.u_dim(:,i) = u2T(DATA_ROMPC.u(:,i), PARAMS);
end
DATA_LQR.u_dim = zeros(size(DATA_LQR.u));
for i = 1:DATA_LQR.T
    DATA_LQR.u_dim(:,i) = u2T(DATA_LQR.u(:,i), PARAMS);
end

%% Plotting state profiles
figure('color',[1,1,1],'Position', [1, 1, 900,1200]); hold on;
% title('Title','Interpreter','latex','Fontsize',22);

[Cstar, Tstar] = x2CT(zeros(nf,1), PARAMS);
pos = linspace(0, 1, PARAMS.P.N);
dt = 0.01;
ks = [1,301, 601,1001];
fsize = 14;
for i = 1:length(ks)
    [C, T] = x2CT(DATA_ROMPC.xf(:,ks(i)), PARAMS);
    title_str = sprintf('$t = %d$ s', (ks(i) - 1)*dt);
    
    subplot(length(ks), 2, 2*(i-1) + 1); hold on;
    plot(pos, T, 'color', colors{2}, 'Linewidth',1);
    plot(pos, Tstar, 'color', colors{2}, 'linestyle','--');
    plot(pos, 395*ones(PARAMS.P.N,1),'k--')
    
    if i == length(ks)
        xlabel('Length, [m]', 'FontSize',fsize,'Interpreter','latex');
    end
    ylabel('Temperature [K]', 'FontSize',fsize,'Interpreter','latex');
    title(title_str,'Interpreter','latex', 'FontSize',fsize);
    
    subplot(length(ks), 2, 2*i); hold on;
    plot(pos, C, 'b', 'Linewidth',1);
    plot(pos, Cstar, 'b', 'linestyle','--');
    
    if i == length(ks)
        xlabel('Length, [m]', 'FontSize',fsize,'Interpreter','latex');
    end
    ylabel('Concentration [mol/L]', 'FontSize',fsize,'Interpreter','latex');
    title(title_str,'Interpreter','latex', 'FontSize',fsize);
end       
% filename = strcat(figpath, 'xf_tubular_reactor');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')

%% Plotting state profiles at specific times, LQR comparison
fsize = 14;
[Cstar, Tstar] = x2CT(zeros(nf,1), PARAMS);
[C0, T0] = x2CT(DATA_ROMPC.xf(:,1), PARAMS);
k = 400;
dt = 0.01;
title_str = sprintf('$t = %d$ s', k*dt);
[Ck, Tk] = x2CT(DATA_ROMPC.xf(:,k), PARAMS);
[Ck_lqr, Tk_lqr] = x2CT(DATA_LQR.xf(:,k), PARAMS);
pos = linspace(0, 1, PARAMS.P.N);

figure('color',[1,1,1],'Position', [1, 1, 900, 300]); hold on;
subplot(1, 2, 1); hold on;
plot(pos, Tk, 'color', colors{2}, 'Linewidth',1,'marker','o','markerindices',[1:30:300]);
plot(pos, Tk_lqr, 'color',colors{3}, 'Linewidth',1,'marker','x','markerindices',[1:30:300]);
plot(pos, 395*ones(PARAMS.P.N,1),'k--')
xlabel('Length, [m]', 'FontSize',fsize,'Interpreter','latex');
ylabel('Temperature [K]', 'FontSize',fsize,'Interpreter','latex');
title(title_str,'Interpreter','latex', 'FontSize',fsize);
legend({'ROMPC', 'ROLQR'}, 'Interpreter','latex',...
            'FontSize',fsize,'Location','best','Orientation','horizontal');
legend('boxoff');

subplot(1, 2, 2); hold on;
plot(pos, Ck, 'color', [0.75, 0, 0.75], 'Linewidth',1,'marker','o','markerindices',[1:30:300]);
plot(pos, Ck_lqr, 'color',[0 0.6 0.3], 'Linewidth',1,'marker','x','markerindices',[1:30:300]);
% plot(pos, Cstar, 'color', [0.75, 0, 0.75], 'linestyle','--');
xlabel('Length, [m]', 'FontSize',fsize,'Interpreter','latex');
ylabel('Concentration [mol/L]', 'FontSize',fsize,'Interpreter','latex');
title(title_str,'Interpreter','latex', 'FontSize',fsize);
legend({'ROMPC', 'ROLQR'}, 'Interpreter','latex',...
            'FontSize',fsize,'Location','best','Orientation','horizontal');
legend('boxoff');
% filename = strcat(figpath, 'xf_tubular_reactor_compare');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')

%% Plot Control
figure('color',[1,1,1],'Position', [1, 1, 900,600]); hold on;
% title('Tubular Reactor','Interpreter','latex','Fontsize',22);
uUB = 395;
uLB = 300;
plot(DATA_ROMPC.u_dim(1,:),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u_dim(2,:),'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1);
plot(DATA_ROMPC.u_dim(3,:),'color',colors{3},'marker',markers{3},'markerindices',mark_idx,'Linewidth',1);
% plot(DATA_LQR.u_dim(1,:),'color',colors{1},'marker',markers{1},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
% plot(DATA_LQR.u_dim(2,:),'color',colors{2},'marker',markers{2},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
% plot(DATA_LQR.u_dim(3,:),'color',colors{3},'marker',markers{3},'markerindices',mark_idx,'Linewidth',1,'Linestyle','--');
% plot(DATA.u(4,:),'color',colors{4},'marker',markers{4},'markerindices',mark_idx,'Linewidth',1);
% plot(DATA.u(5,:),'color',colors{5},'marker',markers{5},'markerindices',mark_idx,'Linewidth',1);
plot(uUB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
plot(uLB*ones(1,DATA_ROMPC.T), 'color','k','Linewidth',1,'Linestyle','--');
ylim([295, 400]);
xlabel('Time, [s]','Interpreter','latex','FontSize',22);
ylabel('Jacket Temperatures, [K]','Interpreter','latex','FontSize',22);
legend({'$T_{J1}$', '$T_{J2}$', '$T_{J3}$', 'constraints'}, 'Interpreter','latex',...
            'FontSize',18,'Location','east','Orientation','vertical');
legend('boxoff');
% filename = strcat(figpath, 'u_tubular_reactor');
% export_fig(filename, '-png', '-m2')
% export_fig(filename, '-pdf')



