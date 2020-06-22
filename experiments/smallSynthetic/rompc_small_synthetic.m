% clear; clc; close all;
path = fileparts(which('rompc_small_synthetic.m'));

% Load system
[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = smallSynthetic(false);
nf = size(FOM.Af, 1);
n = size(ROM.A, 1);
m = size(FOM.Bf, 2);
p = size(FOM.Cf, 1);
o = size(FOM.Hf, 1);
opt.discrete = true;

N = 20;
[ROMPC, XF, SP, Zbar, Ubar, PROB] = buildROMPC(ROM, Z, U, EBOUND, N, opt);

%%
% ROMPC
T = [0, 1];
r = [25];
[xfss, uss, xbarss, ubarss, xhatss] = computeSteadyStateDiscreteTime(FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r);

X0.xf = xfss;
X0.xhat = xhatss;
X0.xbar = xbarss;
T = 50;
opt.NOISE = NOISE; % to include uniform random bounded noise
[DATA_ROMPC] = simulateROMPC(FOM, ROM, CTRL, PROB, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_ROMPC.mat'),'DATA_ROMPC');

% ROLQR
[K,~,~] = dlqr(ROM.A,-ROM.B,ROM.Q,ROM.R);
[DATA_LQR] = simulateROMPC(FOM, ROM, CTRL, K, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_LQR.mat'),'DATA_LQR');


















