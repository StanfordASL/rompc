clear; clc;
path = fileparts(which('mor_sup_diffuser.m'));
[FOM] = supersonicDiffuser();

% Parse out the CFD model and ignore the disturbance part
Af = FOM.Af(1:end-1, 1:end-1);
Bf = FOM.Bf(1:end-1, :);
Bfw = FOM.Af(1:end-1, end);
Cf = FOM.Cf(:,1:end-1);
Hf = FOM.Hf(:,1:end-1);

% Compute ROM
n = 10;
[A, ~, ~, W, V, S] = balancedTruncationDiscrete(Af, Bf, Hf, n);
B = W'*Bf;
Bw = W'*Bfw;
C = Cf*V;
H = Hf*V;

% Now reassemble to include the disturbance terms
A = [A, Bw; zeros(1,n), FOM.Af(end,end)];
B = [B; 0];
Bw = [zeros(n,1); 1];
C = [C, 0];
H = [H, 0; zeros(1,n), 1];
V = [V, zeros(size(V,1), 1); zeros(1,n), 1];
W = [W, zeros(size(W,1), 1); zeros(1,n), 1];

ROM = struct('A', A, 'B', B, 'Bw', Bw, 'C', C, 'H', H, 'V', V, 'W', W, 'S', S);
save(strcat(path, '/data/ROM.mat'), 'ROM');