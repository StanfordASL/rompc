clear; clc;
path = fileparts(which('mor_large_synthetic.m'));
[FOM] = largeSynthetic(true);

% Compute ROM
n = 6;
[A, B, H, W, V, S] = balancedTruncation(FOM.Af, FOM.Bf, FOM.Hf, n, false);
C = FOM.Cf*V;
Bw = [];

ROM = struct('A', A, 'B', B, 'Bw', Bw, 'C', C, 'H', H, 'V', V, 'W', W, 'S', S);
save(strcat(path, '/data/ROM.mat'), 'ROM');