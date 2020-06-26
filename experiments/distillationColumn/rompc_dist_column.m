clear; clc; close all;
path = fileparts(which('rompc_dist_column.m'));

% Load system
[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = distillationColumn(false);
nf = size(FOM.Af, 1);
n = size(ROM.A, 1);
m = size(FOM.Bf, 2);
p = size(FOM.Cf, 1);
o = size(FOM.Hf, 1);
opt.discrete = true;

N = 50;
opt.solver = 'mosek';
[ROMPC, Xf, SP, Zbar, Ubar] = buildROMPC(ROM, Z, U, EBOUND, N, opt);

% Initial conditions (make sure they satisfy initial conditions on sets)
T = [eye(4), zeros(4,4)];
r = [-0.005; -0.009; -0.75; 0.25];
[xfss, uss, xbarss, ubarss, xhatss] = computeSteadyStateDiscreteTime(FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r);

% Simulate ROMPC
X0.xf = xfss;
X0.xhat = xhatss;
X0.xbar = xbarss;
T = 150;
opt.NOISE = NOISE;
[DATA_ROMPC] = simulateROMPC(FOM, ROM, CTRL, ROMPC.prob, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_ROMPC.mat'),'DATA_ROMPC');

% Simulate ROLQR
[K,~,~] = dlqr(ROM.A,-ROM.B,ROM.Q,ROM.R);
[DATA_LQR] = simulateROMPC(FOM, ROM, CTRL, K, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_LQR.mat'),'DATA_LQR');















