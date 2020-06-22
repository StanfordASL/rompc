clear; clc;
path = fileparts(which('mor_small_synthetic.m'));
[FOM] = smallSynthetic(true);

% Compute ROM
n = 4;
[A, B, C, W, V, S] = balancedTruncationDiscrete(FOM.Af, FOM.Bf, FOM.Cf, n);
H = FOM.Hf*V;
Bw = W'*FOM.Bfw;

ROM = struct('A', A, 'B', B, 'Bw', Bw, 'C', C, 'H', H, 'V', V, 'W', W, 'S', S);
save(strcat(path, '/data/ROM.mat'), 'ROM');