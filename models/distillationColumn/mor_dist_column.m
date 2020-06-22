clear; clc;
path = fileparts(which('mor_dist_column.m'));
[FOM, ~, ~, Z, U, NOISE] = distillationColumn(true);

% Add scaling to inputs and outputs
idxz = find(Z.A*ones(8,1) == 1);
Dz = diag(Z.b(idxz)); % max output errors (scalings)
idxu = find(U.A*ones(4,1) == 1);
Du = diag(U.b(idxu)); % max inputs (scalings)
Hfscaled = inv(Dz)*FOM.Hf;
Bfscaled = FOM.Bf*Du;

n = 12;
[A, ~, ~, W, V, S] = balancedTruncationDiscreteUnstable(FOM.Af, Bfscaled, Hfscaled, n);
B = W'*FOM.Bf;
Bw = W'*FOM.Bfw;
C = FOM.Cf*V;
H = FOM.Hf*V;

ROM = struct('A', A, 'B', B, 'Bw', Bw, 'C', C, 'H', H, 'V', V, 'W', W, 'S', S);
save(strcat(path, '/data/ROM.mat'), 'ROM');

