function [constraints, objective, xbar, ubar, var] = setupOptimalControl(A, B, P, Q, R, H, Z, U, Xf, N, xss, uss)
%[constraints, objective, xbar, ubar, var] = setupOptimalControl(A, B, P, Q, R, H, Z, U, Xf, N, xss, uss)
%
%Sets up the optimal control problem at time k with quadratic cost and nominal dynamics
%and returns the constraints list and objective.
%
% Inputs:
%   A: dynamics matrix
%   B: control matrix
%   P: terminal cost matrix
%   Q: state cost matrix (dx'Qdx)
%   R: control cost matrix (du'Rdu)
%   H: state to performance variable matrix (z = Hx)
%   Z: performance variable constraint Polyhedral set
%   U: control constraint Polyhedral set
%   Xf: terminal state constraint Polyhedral set
%   N: time horizon
%   xss: target x value
%   uss: target u value
%
% Outputs:
%   constraints: list of Yalmip constraints
%   objective: Yalmip objective function
%   xbar: Yalmip state variables
%   ubar: Yalmip control variables
%   var: a struct of information

% Here we setup the nominal MPC problem for tracking setpoint
n = size(A, 1); % reduced state dimension
m = size(B, 2); % control dimension

% Define decison variables
xbar = {};
ubar = {};
for i = 1:N
    xbar{i} = sdpvar(n, 1);
    ubar{i} = sdpvar(m, 1);
end
xbar{N+1} = sdpvar(n, 1);

% Quadratic objective function
objective = 0;
for i = 1:N
    objective = objective + (xbar{i} - xss)'*Q*(xbar{i} - xss) + ...
            (ubar{i} - uss)'*R*(ubar{i} - uss);
end
objective = objective + (xbar{N+1} - xss)'*P*(xbar{N+1} - xss); % terminal cost

% Constraints
constraints = [];

% Add dynamics constraint
for i = 1:N
    constraints = [constraints, xbar{i+1} == A*xbar{i} + B*ubar{i}];
end

% Add tightened state and control constraints
for i = 1:N
    constraints = [constraints, Z.A*H*xbar{i} <= Z.b];
    constraints = [constraints, U.A*ubar{i} <= U.b];
end

% Add terminal constraint
if isa(Xf, 'Polyhedron')
    constraints = [constraints, Xf.A*xbar{N+1} <= Xf.b];
end

% Log some info about the problem
var.N = N;
var.n = n;
var.m = m;

end

