clear; clc;
path = fileparts(which('control_small_synthetic.m'));

[FOM, ROM, ~, Z, U, NOISE] = smallSynthetic(false);
opts.discrete = true;
opts.with_noise = false;

% Compute controller
opts.method = 'ric';
opts.Wz = diag([1,10]);
[CTRL.K, CTRL.L] = computeControllerGains(FOM, ROM, opts);
save(strcat(path, '/data/CTRL.mat'), 'CTRL');

[ERROR] = buildErrorSystem(FOM, ROM, CTRL, Z, U, NOISE.W, NOISE.V);
max(abs(eig(ERROR.Ae)))



