clear; clc;
path = fileparts(which('mor_heatflow.m'));
[FOM] = heatflow(true);

% Compute ROM
n = 21;
[A, B, C, W, V, S] = balancedTruncationDiscreteUnstable(FOM.Af, FOM.Bf, FOM.Cf, n);
Bw = W'*FOM.Bfw;
H = FOM.Hf*V;

ROM = struct('A', A, 'B', B, 'Bw', Bw, 'C', C, 'H', H, 'V', V, 'W', W, 'S', S);
save(strcat(path, '/data/ROM.mat'), 'ROM');