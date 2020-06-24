function [x, u, diagnostics] = solveOptimalControl(PROB, x0, opt)
%[x, u, diagnostics] = solveOptimalControl(PROB, x0, opt)
%
%Solves the optimal control problem at time k given x0 and returns the
%optimal x_k and u_k sequences.
% Inputs:
%   PROB: struct containing
%       - constraints: list of Yalmip constraints for the optimization problem
%       - objective: Yalmip objective for the optimization problem
%       - x: Yalmip state variables
%       - u: Yalmip control variables
%       - var: Yalmip stuff
%       - solver: solver to use
%        This struct is output by the function buildMPC() automatically
%   x0: initial condition for the MPC
%   opt: struct to add additional options
%       - verbose: true or false
%
% Outputs:
%   x: optimal state sequence at time k
%   u: optimal control sequence at time k
%   diagnostics: solver output

if ~isfield(opt, 'verbose')
    opt.verbose = false;
end

% Add constraint for initial condition
constraints = [PROB.constraints, PROB.x{1} == x0];

% Solve the problem
opt = sdpsettings('verbose', opt.verbose, 'solver', PROB.solver, 'savesolveroutput', 1);
diagnostics = optimize(constraints, PROB.objective, opt);

x = zeros(PROB.var.n, PROB.var.N+1);
u = zeros(PROB.var.m, PROB.var.N);
for i = 1:PROB.var.N
    x(:,i) = value(PROB.x{i});
    u(:,i) = value(PROB.u{i});
end
x(:, PROB.var.N+1) = value(PROB.x{PROB.var.N+1});

end

