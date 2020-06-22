clear; clc;
path = fileparts(which('control_sup_diffuser.m'));
[FOM, ROM, ~, Z, U, NOISE] = supersonicDiffuser();
opts.discrete = true;

% Compute controller
opts.method = 'ric';
opts.Wz = diag([1,0]);
opts.Wu = [0.5];
[CTRL.K, CTRL.L] = computeControllerGains(FOM, ROM, opts);
save(strcat(path, '/data/CTRL.mat'), 'CTRL');



