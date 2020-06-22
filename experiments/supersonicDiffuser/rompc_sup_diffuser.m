%clear; clc; close all;
path = fileparts(which('rompc_sup_diffuser.m'));

% Load system
%[FOM, ROM, CTRL, Z, U, NOISE, EBOUND, params] = supersonicDiffuser();
nf = size(FOM.Af, 1);
n = size(ROM.A, 1);
m = size(FOM.Bf, 2);
p = size(FOM.Cf, 1);
o = size(FOM.Hf, 1);
opt.discrete = true;

N = 20;
EBOUND.Dz_2(3:4) = 0; % Don't need actually need constraints on the estimated d
[ROMPC, XF, SP, Zbar, Ubar, PROB] = buildROMPC(ROM, Z, U, EBOUND, N, opt);

% ROMPC
T = [1, 0];
r = [0.05];
[xfss, uss, xbarss, ubarss, xhatss] = computeSteadyStateDiscreteTime(FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r);

X0.xf = xfss;
X0.xhat = xhatss;
X0.xbar = xbarss;
T = 1000;

% Create Gaussian noise
t = linspace(0, params.dt*T, T);
t0 = t(end)/2;
wf_profile = params.M*exp(-params.alpha*(t-t0).^2);
opt.NOISE.wf = [wf_profile(2:end) - 0.95*wf_profile(1:end-1), 0];
opt.NOISE.V = NOISE.V; % to include uniform random bounded noise

% Simulate with the Gaussian noise profile
[DATA_ROMPC] = simulateROMPC(FOM, ROM, CTRL, ROMPC, Zbar, T, X0, opt);
DATA_ROMPC.wf = opt.NOISE.wf;
DATA_ROMPC.dt = params.dt;
save(strcat(path, '/data/DATA_ROMPC.mat'),'DATA_ROMPC');


% Simulate with the Gaussian noise profile, but an LQR controller that does
% not include any disturbance estimation
ROM2.A = ROM.A(1:end-1,1:end-1);
ROM2.B = ROM.B(1:end-1,:);
ROM2.Bw = ROM.Bw(1:end-1,:);
ROM2.C = ROM.C(:,1:end-1);
ROM2.H = ROM.H(1,1:end-1);
ROM2.Q = ROM.Q(1:end-1,1:end-1);
ROM2.R = ROM.R;
ctrl_opts.discrete = 'true';
ctrl_opts.method = 'ric';
ctrl_opts.Wu = [0.5];
[CTRL2.K, CTRL2.L] = h2optController(FOM, ROM2, ctrl_opts);
[K,~,~] = dlqr(ROM2.A,-ROM2.B,ROM2.Q,ROM2.R);

% Simulate LQR case
X02.xf = xfss;
X02.xhat = xhatss(1:end-1);
X02.xbar = xbarss(1:end-1);
[DATA_LQR] = simulateROMPC(FOM, ROM2, CTRL2, K, Zbar, T, X02, opt);
DATA_LQR.wf = opt.NOISE.wf;
DATA_LQR.dt = params.dt;
save(strcat(path, '/data/DATA_LQR.mat'),'DATA_LQR');



















