clear; clc;
path = fileparts(which('control_tubular_reactor.m'));

[FOM, ROM, ~, Z, U, NOISE] = tubularReactor(false);
opts.continuous = true;

% Compute controller
opts.method = 'ric';
[CTRL.K, CTRL.L] = computeControllerGains(FOM, ROM, opts);
save(strcat(path, '/data/CTRL.mat'), 'CTRL');

[ERROR] = buildErrorSystem(FOM, ROM, CTRL, Z, U, NOISE.W, NOISE.V);
max(real(eig(ERROR.Ae)))

