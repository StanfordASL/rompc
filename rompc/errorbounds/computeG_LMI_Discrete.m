function [G, diagnostics] = computeG_LMI_Discrete(P, opt)
%[G, diagnostics] = computeG_LMI_Discrete(P, opt)
%
%This function seeks to compute a positive definite matrix G such that the
%quantity gamma is minimized, where ||P||_G <= gamma.
%
% Inputs:
%   P: n x n matrix that is Schur stable
%   opt: opt.solver specifies the solver to use
%
% Returns:
%   G: n x n diagonal, positive definite matrix

if ~isfield(opt, 'solver')
    opt.solver = 'mosek';
end

n = size(P, 1); % dimension of G

% Formulate optimization problem
G = sdpvar(n, n);
gamma = sdpvar(1);
ops = sdpsettings('verbose', 1, 'solver', opt.solver);

% Solve via bisection search
diagnostics = bisection([P'*G*P <= gamma^2*G, G >= 1, gamma >= 0], gamma, ops);
if diagnostics.problem ~= 0
    disp('Solver failed.');
    disp(diagnostics.info)
end
G = value(G);

end

