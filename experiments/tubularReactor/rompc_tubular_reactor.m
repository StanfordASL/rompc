clear; clc; close all;
path = fileparts(which('rompc_tubular_reactor.m'));

% Load system
[FOM, ROM, CTRL, Z, U, NOISE, EBOUND, PARAMS] = tubularReactor(false);
opt.continuous = true;

N = 50;
opt.xf_datapath = strcat(path, '/data/XF.mat');
opt.solver = 'mosek';
[ROMPC, Xf, SP, Zbar, Ubar] = buildROMPC(ROM, Z, U, EBOUND, N, opt);

% Initial conditions
TJ1 = 300;
TJ2 = 300;
TJ3 = 300;
[C0, T0] = reactorSteadyState(TJ1, TJ2, TJ3, PARAMS.P);
[X0.xf] = CT2x(C0, T0, PARAMS);
X0.xbar = ROM.W'*X0.xf;
X0.xhat = X0.xbar;

% Simulate
T = 15; % seconds
opt.sim_dt = 0.01;
opt.rompc_dt = ROMPC.dt;
opt.save_xf = true;
[DATA_ROMPC] = simulateROMPC(FOM, ROM, CTRL, ROMPC.rompc, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_ROMPC.mat'),'DATA_ROMPC');

% ROLQR
[K,~,~] = lqr(ROM.A,-ROM.B,ROM.Q,ROM.R);
[DATA_LQR] = simulateROMPC(FOM, ROM, CTRL, K, Zbar, T, X0, opt);
save(strcat(path, '/data/DATA_LQR.mat'),'DATA_LQR');


















