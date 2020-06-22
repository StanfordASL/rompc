function [MPC, PROB] = buildMPC(A, B, P, Q, R, H, Z, U, Xf, N, xss, uss, opt)
%[MPC, PROB] = buildMPC(A, B, P, Q, R, H, Z, U, Xf, N, xss, uss, opt)
%
%Defines the MPC optimal control problem. The output is an object with a pre-compiled
%lower-level numerical format that reduces MPC overhead. For MPC, the
%resulting control can be computed via
%
% u(k) = MPC{x(k)};
%
% Inputs:
%   A: dynamics matrix
%   B: control matrix
%   P: terminal cost matrix
%   Q: state cost matrix
%   R: control cost matrix
%   H: state to performance variable matrix
%   Z: performance variable constraint Polyhedral set
%   U: control constraint Polyhedral set
%   Xf: terminal state constraint (Polyhedral set)
%   N: time horizon
%   xss: target x value
%   uss: target u value
%   opt: struct to add additional options
%       - solver: specify solver for YALMIP (string)
%
% Outputs:
%   ROMPC: Yalmip optimizer (compiled low level object for solving problem)
%   PROB: optimization objects (Yalmip) to call optimize

% Check positive definiteness of Q and R
mineig = min(real(eig(Q)));
if mineig <= 1e-6 && mineig >= 0
    Q = Q + 1e-6*eye(size(Q,1));
elseif mineig < 0
    error('Q needs to be positive definite.');
end

mineig = min(real(eig(R)));
if mineig <= 1e-6 && mineig >= 0
    R = R + 1e-6*eye(size(Q,1));
elseif mineig < 0
    error('R needs to be positive definite.');
end

% Set up the optimal control problem
[constraints, objective, x, u, var] = setupOptimalControl(A, B, P, Q, R, H, Z, U, Xf, N, xss, uss);

% Set up the solver
if isfield(opt, 'solver')
    solver = opt.solver;
else
    solver = 'mosek';
end
solver_opt = sdpsettings('verbose', 0, 'solver', solver);

% optimizer(constraints, objective, options, parameters, wantedvariables)
MPC = optimizer(constraints, objective, solver_opt, {x{1}}, {u{1}, x{2}});

PROB.constraints = constraints;
PROB.objective = objective;
PROB.x = x;
PROB.u = u;
PROB.var = var;
PROB.solver = solver;

end

