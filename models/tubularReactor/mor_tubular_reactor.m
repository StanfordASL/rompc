clear; clc;
path = fileparts(which('mor_tubular_reactor.m'));
[FOM] = tubularReactor(true);

% Compute ROM
n = 30;
[A, B, ~, W, V, S] = balancedTruncationContinuous(FOM.Af, FOM.Bf, [FOM.Cf; FOM.Hf], n);
C = FOM.Cf*V;
H = FOM.Hf*V;
Bw = W'*FOM.Bfw;

ROM = struct('A', A, 'B', B, 'Bw', Bw, 'C', C, 'H', H, 'V', V, 'W', W, 'S', S);
save(strcat(path, '/data/ROM.mat'), 'ROM');