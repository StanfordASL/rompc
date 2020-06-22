clear; clc;
path = fileparts(which('control_dist_column.m'));

[FOM, ROM, ~, Z, U, NOISE] = distillationColumn(false);
opts.discrete = true;
opts.with_noise = false;

% Add scaling to inputs and outputs
idxz = find(Z.A*ones(8,1) == 1);
Dz = diag(Z.b(idxz)); % max output errors (scalings)
idxu = find(U.A*ones(4,1) == 1);
Du = diag(U.b(idxu)); % max inputs (scalings)
opts.Wz = inv(Dz);
opts.Wu = Du;

% Compute controller
opts.method = 'ric';
[CTRL.K, CTRL.L] = computeControllerGains(FOM, ROM, opts);
save(strcat(path, '/data/CTRL.mat'), 'CTRL');

% Check largest eigenvalue
[ERROR] = buildErrorSystem(FOM, ROM, CTRL, Z, U, NOISE.W, NOISE.V);
max(abs(eig(ERROR.Ae)))


