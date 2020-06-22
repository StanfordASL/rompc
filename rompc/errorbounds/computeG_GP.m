function [G, solver] = computeG_GP(X, Y)
%[G, solver] = computeG_GP(X, Y)
%
%This function seeks to minimize the quantitiy \prod_i ||G^1/2*X_i||_2 * ||Y_i*G^-1/2||_2
%by noting that ||A||_2 <= ||A||_F and therefore minimizing the square of
%the upper bound: \prod_i (||G^(1/2)*X_i||_F * ||Y_i*G^(-1/2)||_F)^2. This is solved
%via geometric programming by assuming G is diagonal and PD.
%
% Inputs:
%   X: cell array of X_i matrices, which are n x m_xi in size
%   Y: cell array of X_i matrices, which are n_yi x n in size
%
% Returns:
%   G: n x n diagonal, positive definite matrix

n = size(X{1}, 1); % dimension of G
num_sets = length(X);

% Formulate geometric program.
J = 1;
g = sdpvar(n,1);
fprintf('Building geometric program.\n');
for i = 1:num_sets
    ai = zeros(n, 1);
    bi = zeros(n, 1);
    
    % ai(j) is the sum of the squared values of the jth row of Xi
    % bi(j) is the sum of the squared values of the jth column of Yi
    for j = 1:n
        ai(j) = X{i}(j,:)*X{i}(j,:)';
        bi(j) = Y{i}(:,j)'*Y{i}(:,j);
    end
    
    % Add to cost function
    J = J*(ai'*g)*(bi'*g.^(-1));
end

% Solve geometric program.
fprintf('Optimizing G via geometric programming.\n');
ops = sdpsettings('verbose', 0, 'solver', 'mosek');
solver = optimize([g >= 0], J, ops);
if solver.problem ~= 0
    fprintf('Geometric program failed to find solution in optimizeG.\n');
    G = zeros(n);
else
    G = diag(value(g));
end
end

