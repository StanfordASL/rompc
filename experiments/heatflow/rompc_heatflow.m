clear; clc; close all;
path = fileparts(which('rompc_heatflow.m'));

% Load system
[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = heatflow(false);
nf = size(FOM.Af,1);
n = size(ROM.A,1);
opt.discrete = true;

N = 30;
opt.setpoint.T = eye(size(ROM.H,1));
opt.setpoint.r = [0; 1; -0.1; 0.1];
opt.setpoint.FOM = FOM;
opt.setpoint.CTRL = CTRL;
opt.solver = 'cplex';
[ROMPC, Xf, SP, Zbar, Ubar] = buildROMPC(ROM, Z, U, EBOUND, N, opt);

% ROMPC
xfss = zeros(nf,1);
xbarss = zeros(n,1);
xhatss = zeros(n,1);

X0.xf = xfss;
X0.xhat = xhatss;
X0.xbar = xbarss;
T = 500;
[DATA_ROMPC] = simulateROMPC(FOM, ROM, CTRL, ROMPC.rompc, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_ROMPC.mat'),'DATA_ROMPC');



















