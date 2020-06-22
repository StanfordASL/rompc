function [EexpAtB_approx] = recursiveImplicitEuler(A, B, E, dt, N)
%[EexpAtB_approx] = recursiveImplicitEuler(A, B, E, dt, N)
%
% For the continuous time system xd = A*x + B*u this computes approximations 
% to the set of matrices [E*B, E*e^(Adt)*B, ... , E*e^(Adt*N)*B].
% Rather than computing the matrix exponential, this approach performs an
% implicit Euler integration of the continuous time dynamics.
%
% Inputs:
%   A: continuous time system matrix
%   B: continuous time system matrix
%   E (optional): output matrix or [] to use I
%   dt: time discretization
%   N: number of time steps for recursion
%
% Outputs:
%   EexpAtB_approx: cell array [E*B, E*e^(Adt)*B, ... , E*e^(Adt*N)*B] of
%                  approximate matrices

fprintf('Using the implicit Euler method.\n');
EexpAtB_approx = cell(1, N+1); % cell array to store matrices
n = size(A, 1);
m = size(B, 2);

if isempty(E)
    E = eye(n);
end

% Compute initial condition for recursion
V = B;
EexpAtB_approx{1} = E*V;

% Computing preconditioner
fprintf('Computing preconditioner.\n');
S = speye(n) - dt*sparse(A);
setup.type = 'nofill';
[M1, M2] = ilu(S, setup);

% Perform recursion
fprintf('Precomputing some matrices.      ');
maxiter = min(1000, n);
for i = 2:N+1
    V_prev = V;
    for j = 1:m
        [V(:,j), ~, ~] = gmres(S, V_prev(:,j), [], [], maxiter, M1, M2, []);
    end
    EexpAtB_approx{i} = E*V;
    fprintf('\b\b\b\b\b%3.0f%% ', round(100*(i/(N+1))));
end
fprintf('\n');

end

