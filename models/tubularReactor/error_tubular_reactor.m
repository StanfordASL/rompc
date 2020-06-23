clear; clc;
path = fileparts(which('error_tubular_reactor.m'));

[FOM, ROM, CTRL, Z, U, NOISE] = tubularReactor(false);
opts.continuous = true;

opt.G_method = 'Lyap';
opts.dt = 0.05; % for continuous time
opt.mat_exp_method = 'expm';
opt.Cr_method = 'upperbound';
opt.Cw_method = 'upperbound';
tau = 60;
opts.normbound_datapath = strcat(path, '/data/NORMBOUND.mat');
opts.errormats_datapath = strcat(path, '/data/ERRORMATS.mat');
opts.xbar_datapath = strcat(path, '/data/XBAR.mat');
[EBOUND] = computeErrorBounds(FOM, ROM, CTRL, Z, U, NOISE.W, NOISE.V, tau, opts);
save(strcat(path, '/data/EBOUND.mat'), 'EBOUND');

