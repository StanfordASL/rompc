clear; clc;
path = fileparts(which('control_heatflow.m'));

[FOM, ROM, ~, Z, U, NOISE] = heatflow(false);
o = size(FOM.Hf,1);
opts.discrete = true;

% Compute controller
opts.method = 'ric';
opts.Wz = 100*eye(o);
[CTRL.K, CTRL.L] = computeControllerGains(FOM, ROM, opts);
save(strcat(path, '/data/CTRL.mat'), 'CTRL');

[ERROR] = buildErrorSystem(FOM, ROM, CTRL, Z, U, NOISE.W, NOISE.V);
max(abs(eig(ERROR.Ae)))



