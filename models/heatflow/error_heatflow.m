clear; clc;
path = fileparts(which('error_heatflow.m'));

[FOM, ROM, CTRL, Z, U, NOISE] = heatflow(false);
opts.discrete = true;

opts.G_method = 'BatchGP';
opts.Cr_method = 'hyperrect';
opts.Cw_method = 'hyperrect';
tau = 2000;
opts.Xbar_tau = 500;
opts.normbound_datapath = strcat(path, '/data/NORMBOUND.mat');
opts.errormats_datapath = strcat(path, '/data/ERRORMATS.mat');
opts.xbar_datapath = strcat(path, '/data/XBAR.mat');
[EBOUND] = computeErrorBounds(FOM, ROM, CTRL, Z, U, NOISE.W, NOISE.V, tau, opts);
save(strcat(path, '/data/EBOUND.mat'), 'EBOUND');

