function [proj] = correctionProjection(H, Zbar, opt)
%[proj] = correctionProjection(H, Zbar, opt)
%
%YALMIP when outputing the solution sometimes run into numerical issues
%where the next iterate is not initially feasible (this occurs when
%solution rides a constraint boundary). So attempt to correct for this my
%reprojecting the solution into the interior of Zbar.

if ~isfield(opt, 'solver')
    opt.solver = 'cplex';
end

n = size(H,2);
x = sdpvar(n,1);
xcur = sdpvar(n,1);
constraints = [Zbar.A*H*x <= Zbar.b - .0001];
objective = norm(x - xcur, 1);
ops = sdpsettings('verbose', 0, 'solver', opt.solver, 'savesolveroutput', 1);
proj = optimizer(constraints, objective, ops, xcur, x);
end

