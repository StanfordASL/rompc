function [G, alpha, diagnostics] = computeG_LMI_Continuous(P, opt)
%[G, alpha, diagnostics] = computeG_LMI_Continuous(P, opt)
%
%This function seeks to compute a positive definite matrix G such that the
%quantity alpha is minimized, where ||e^Pt||_G <= e^(alpha*t).
%
% Inputs:
%   P: n x n matrix that is Hurwitz
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
alpha = sdpvar(1);
ops = sdpsettings('verbose', 1, 'solver', opt.solver);

% Solve via bisection search
diagnostics = bisection([P'*G + G*P <= 2*alpha*G, G >= 1], alpha, ops);
if diagnostics.problem ~= 0
    disp('Solver failed.');
    disp(diagnostics.info)
end
G = value(G);
alpha = value(alpha);
end

