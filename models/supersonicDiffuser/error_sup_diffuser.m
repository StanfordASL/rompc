clear; clc;
path = fileparts(which('error_sup_diffuser.m'));
[FOM, ROM, CTRL, Z, U, NOISE] = supersonicDiffuser();
opts.discrete = true;

opts.G_method = 'Lyap';
opts.Cr_method = 'hyperrect';
opts.Cw_method = 'hyperrect';
tau = 2000;
opts.normbound_datapath = strcat(path, '/data/NORMBOUND.mat');
opts.errormats_datapath = strcat(path, '/data/ERRORMATS.mat');
opts.xbar_datapath = strcat(path, '/data/XBAR.mat');
[EBOUND] = computeErrorBounds(FOM, ROM, CTRL, Z, U, NOISE.W, NOISE.V, tau, opts);
save(strcat(path, '/data/EBOUND.mat'), 'EBOUND');

% Also compute some without noise
[EBOUND] = computeErrorBounds(FOM, ROM, CTRL, Z, U, 0, 0, tau, opts);
save(strcat(path, '/data/EBOUND_no_noise.mat'), 'EBOUND');



