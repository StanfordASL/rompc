clear; clc;
path = fileparts(which('control_large_synthetic.m'));

[FOM, ROM, ~, Z, U, NOISE] = largeSynthetic(false);
opts.discrete = true;

% Compute controller
opts.method = 'ric';
[CTRL.K, CTRL.L] = computeControllerGains(FOM, ROM, opts);
save(strcat(path, '/data/CTRL.mat'), 'CTRL');

[ERROR] = buildErrorSystem(FOM, ROM, CTRL, Z, U, NOISE.W, NOISE.V);
max(abs(eig(ERROR.Ae)))



