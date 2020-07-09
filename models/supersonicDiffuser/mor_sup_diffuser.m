clear; clc;
path = fileparts(which('mor_sup_diffuser.m'));
[FOM] = supersonicDiffuser();

% Compute ROM
n = 10;
[A, ~, ~, W, V, S] = balancedTruncation(FOM.Af, FOM.Bf, FOM.Hf, n, false);
B = W'*Bf;
Bw = W'*Bfw;
C = Cf*V;
H = Hf*V;

ROM = struct('A', A, 'B', B, 'Bw', Bw, 'C', C, 'H', H, 'V', V, 'W', W, 'S', S);
save(strcat(path, '/data/ROM.mat'), 'ROM');