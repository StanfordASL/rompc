function [prob] = supportFunction(N, M, opt)
%[prob] = supportFunction(N, M, opt)
%
% Builds the support function optimization of a polytope hw(a)
%   max. a'*w   s.t. F*w <= g
%
% Returns the optimizer object, which can be called as
%   [w,flag,~,~,~,diagnostics] = prob(a, F, g);
%
% Inputs:
%   N: size of the vector w and a
%   M: size of the number of constraints
%   opt: use opt.solver to choose solver
%
% Outputs:
%   Jstar: solution to the LP (optimizer = false)
%   diagnostics: solver info (optimizer = false)
%   prob: optimizer object,  (optimizer = true)

opt.solver = 'cplex';

a = sdpvar(N, 1);
w = sdpvar(N, 1);
F = sdpvar(M, N);
g = sdpvar(M, 1);
constraints = [F*w <= g];
objective = -a'*w;
ops = sdpsettings('verbose', 0, 'solver', opt.solver, 'savesolveroutput', 1);
prob = optimizer(constraints, objective, ops, {a, F, g}, w);

end

