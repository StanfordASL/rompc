function [EexpAtB] = recursiveMatExpAction(A, B, E, dt, N)
%[EexpAtB] = recursiveMatExpAction(A, B, E, dt, N)
%
% For the continuous time system xd = A*x + B*u this computes the
% set of matrices [E*B, E*e^(Adt)*B, ... , E*e^(Adt*N)*B].
% Rather than computing the matrix exponential, this function only computes
% actions of the matrix exponential (i.e. e^(Adt)v).
%
% Inputs:
%   A: continuous time system matrix
%   B: continuous time system matrix
%   E (optional): output matrix or [] to use I
%   dt: time discretization
%   N: number of time steps for recursion
%
% Outputs:
%   EexpAtB: cell array [E*B, E*e^(Adt)*B, ... , E*e^(Adt*N)*B]

fprintf('Using the matrix exponential action method.\n');
EexpAtB = cell(1, N+1); % cell array to store matrices
n = size(A, 1);
m = size(B, 2);
M_A = [A, B; sparse(m, n + m)];

if isempty(E)
    E = speye(n);
end

fprintf('Precomputing some matrices.      ');
% Compute initial condition for recursion
V = B;
EexpAtB{1} = E*V;

% Perform recursion
for i = 2:N+1
    V_prev = V;
    [V,~,~,~,~] = expmv(dt, A, V_prev);
    EexpAtB{i} = E*V;
    fprintf('\b\b\b\b\b%3.0f%% ', round(100*(i/(N+1))));
end
fprintf('\n');

end

