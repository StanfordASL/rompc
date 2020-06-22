 clear; clc; close all;

[FOM, ROM, ~, Z, U, Wnoise, Vnoise] = randomSmallSystem(false);
nf = size(FOM.Af,1);
n = size(ROM.A,1);
m = size(FOM.Bf,2);
p = size(FOM.Cf,1);
o = size(FOM.Hf,1);
FOM.Bfw = [];

opts.discrete = true;

% Compute LQR controller
opts.method = 'ric';
[K, L] = h2optController(FOM, ROM, opts);
CTRL_RIC.K = K;
CTRL_RIC.L = L;
[n_ric] = computeErrorH2Norm(FOM, ROM, CTRL_RIC, opts);

% Compute H2SOFO controller
opts.method = 'h2sofo';
opts.warmstart = true;
FOM.Bfw = [];
[K, L] = h2optController(FOM, ROM, opts);
CTRL_H2SOFO.K = K;
CTRL_H2SOFO.L = L;
[n_H2SOFO] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO, opts);

% Compute new controller
opts.method = 'lmi';
[K, L] = h2optController(FOM, ROM, opts);
CTRL_LMI.K = K;
CTRL_LMI.L = L;
[n_LMI] = computeErrorH2Norm(FOM, ROM, CTRL_LMI, opts);

fprintf('LQR: %.4f, H2SOFO: %.4f, LMI: %.4f\n', n_ric, n_H2SOFO, n_LMI);

%% Compute error bounds
tau = 100;
opt = {};
[Xbar] = computeXbar(ROM.A, ROM.B, ROM.H, Z, U, tau, 1, opt);
[Dz_LQR, Du_LQR] =  computeInputEffectBoundDiscrete(ROM, ERROR_LQR, Z, U, Xbar, 0, 0, tau, opt);
[Dz_H2SOFO, Du_H2SOFO] =  computeInputEffectBoundDiscrete(ROM, ERROR_H2SOFO, Z, U, Xbar, 0, 0, tau, opt);
[Dz_LMI, Du_LMI] =  computeInputEffectBoundDiscrete(ROM, ERROR_LMI, Z, U, Xbar, 0, 0, tau, opt);



