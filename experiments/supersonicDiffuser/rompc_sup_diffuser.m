clear; clc; close all;
path = fileparts(which('rompc_sup_diffuser.m'));

% Load system
[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = supersonicDiffuser();
opt.discrete = true;

N = 20;
[ROMPC, Xf, SP, Zbar, Ubar] = buildROMPC(ROM, Z, U, EBOUND, N, opt);

% ROMPC
T = [1];
r = [0.05];
[xfss, uss, xbarss, ubarss, xhatss] = computeSteadyStateDiscreteTime(FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r);

X0.xf = xfss;
X0.xhat = xhatss;
X0.xbar = xbarss;
T = 80;

% Simulate with noise
opt.NOISE = NOISE;
[DATA_ROMPC] = simulateROMPC(FOM, ROM, CTRL, ROMPC.rompc, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_ROMPC.mat'),'DATA_ROMPC');


% Simulate LQR case
[K,~,~] = dlqr(ROM.A, -ROM.B, ROM.Q, ROM.R);
[DATA_LQR] = simulateROMPC(FOM, ROM, CTRL, K, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_LQR.mat'),'DATA_LQR');



















