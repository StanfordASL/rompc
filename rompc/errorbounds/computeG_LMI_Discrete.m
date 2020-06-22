function [G, diagnostics] = computeG_LMI_Discrete(P)
%[G, diagnostics] = computeG_LMI_Discrete(P)
%
%This function seeks to compute a positive definite matrix G such that the
%quantity gamma is minimized, where ||P||_G <= gamma.
%
% Inputs:
%   P: n x n matrix that is Schur stable
%
% Returns:
%   G: n x n diagonal, positive definite matrix

n = size(P, 1); % dimension of G

% Formulate optimization problem
G = sdpvar(n, n);
gamma = sdpvar(1);
ops = sdpsettings('verbose', 1, 'solver', 'mosek');

% Solve via bisection search
diagnostics = bisection([P'*G*P <= gamma^2*G, G >= 1, gamma >= 0], gamma, ops);
if diagnostics.problem ~= 0
    disp('Solver failed.');
    disp(diagnostics.info)
end
G = value(G);

end

