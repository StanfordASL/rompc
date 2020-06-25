clear; clc;
path = fileparts(which('error_aircraft.m'));

[FOM, ROM, CTRL, Z, U, NOISE] = aircraft();
opts.continuous = true;

opts.no_FOM = true;
opts.D2_only = true;
opts.dt = 0.05;
tau = 75;
opts.Xbar_tau = 10;
opts.Xbar_t = opts.Xbar_tau;
opts.errormats_datapath = strcat(path, '/data/ERRORMATS.mat');
opts.xbar_datapath = strcat(path, '/data/XBAR.mat');
[EBOUND] = computeErrorBounds(FOM, ROM, CTRL, Z, U, 0, 0, tau, opts);
save(strcat(path, '/data/EBOUND.mat'), 'EBOUND');


