function [x] = simulateViaBackEuler(A, B, x0, dt, u, w, M1, M2)
%[x] = simulateViaBackEuler(A, B, x0, dt, u, w, M1, M2)
%
% Integrates the continuous time ODE xdot = Ax + Bu + w over N steps where
% N is the number of inputs given in u and w
%
% Inputs:
%   A: n x n matrix
%   B: n x m matrix
%   x0: n x 1 vector, initial condition
%   dt: time step
%   u: m x N matrix of control inputs
%   w: n x N matrix of disturbances
%   M1: preconditioners (optional)
%   M2: preconditioners (optional)
%
% Outputs:
%   x: solution

N = size(u,2); % number of time steps to integrate over
n = size(A,1); % size of the system
x = zeros(n, N+1);
x(:,1) = x0; % initial condition

% Integrate
if issparse(A)
    S = speye(n) - dt*A;
    if ~exist('M1','var') || ~exist('M2','var')
        setup.type = 'nofill';
        [M1, M2] = ilu(S, setup);
    end
else
    S = eye(n) - dt*A;
end
for k = 1:N
    % Solve the implicit Backward Euler step
    if issparse(A)
        [x(:,k+1), flag, res] = gmres(S, x(:,k) + dt*(B*u(:,k) + w(:,k)), [], [], 1000, M1, M2, x(:,k));
        if flag ~= 0
            fprintf('gmres failed to converge to tolerance within max iterations, res=%f\n',res);
        end
    else
        x(:,k+1) = linsolve(S, x(:,k) + dt*(B*u(:,k) + w(:,k)));
    end
end

end

