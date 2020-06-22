function [G] = computeG_BatchGP(X, Y, N, method, maxiter, epsilon, G0)
%[G] = computeG_BatchGP(X, Y, N, method, maxiter, epsilon, G0)
%
%This function seeks to minimize the quantitiy \prod_i ||G^1/2*X_i||_2 * ||Y_i*G^-1/2||_2
%by noting that ||A||_2 <= ||A||_F and therefore minimizing the square of
%the upper bound: \prod_i (||G^(1/2)*X_i||_F * ||Y_i*G^(-1/2)||_F)^2. This is solved
%via geometric programming by assuming G is diagonal and PD.
%
% Inputs:
%   X: cell array of X_i matrices, which are n x m_xi in size
%   Y: cell array of X_i matrices, which are n_yi x n in size
%   N: batch size
%   method: string, approach for creating batch ('sequential' or 'random')
%   maxiter: max number of iterations
%   epsilon: percent improvement termination condition (%)
%   G0: vector of strictly positive values for initial guess
%
% Returns:
%   G: n x n diagonal, positive definite matrix

n = size(X{1}, 1); % dimension of G
num_sets = length(X);

if strcmp(method, 'sequential')
    S_k = [0];
end

% Precompute some stuff
a = zeros(n, num_sets); 
b = zeros(n, num_sets);
for j = 1:num_sets
    % a(i,j) is the sum of the squared values of the ith row of Xj
    % b(i,j) is the sum of the squared values of the ith column of Yj
    for i = 1:n
        a(i,j) = X{j}(i,:)*X{j}(i,:)';
        b(i,j) = Y{j}(:,i)'*Y{j}(:,i);
    end
end

G_k = G0;
g = sdpvar(N,1);
ops = sdpsettings('verbose', 0, 'solver', 'mosek');
iter = 0;
J_old = inf;
while true
    fprintf('i = %d) ', iter + 1);
    
    % Define batch
    if strcmp(method, 'sequential')
        if S_k(end) == n % start over
            S_k = [1:N];
        else % move the window
            first = S_k(end) + 1;
            last = (first + N) - 1;
            if last > n
                S_k = [first:n];
            else
                S_k = [first:last];
            end
        end        
    elseif strcmp(method, 'random')
        num = 0;
        S_k = [];
        while num < N
            cand = randi(n, 1);
            if ~ismember(cand, S_k)
                S_k = [S_k, cand];
                num = num + 1;
            end
        end
    else
        return;
    end
    
    % Formulate geometric program.
    fprintf('Building. ');
    J = 1;
    for j = 1:num_sets
        % Cost function definition
        J1 = 0; J2 = 0;
        for i = 1:n
            if ismember(i, S_k)
                index = find(S_k == i);
                J1 = J1 + g(index)*a(i,j);
                J2 = J2 + (1/g(index))*b(i,j);
            else
                J1 = J1 + G_k(i)*a(i,j);
                J2 = J2 + (1/G_k(i))*b(i,j);
            end
        end
        J = J*J1*J2;
    end

    % Solve geometric program.
    fprintf('Solving. ');
    solver = optimize([g >= 0], J, ops);
    if solver.problem ~= 0
        fprintf('Geometric program failed to find solution in optimizeG.\n');
        G = zeros(n);
    else
        for i=1:length(S_k)
            G_k(S_k(i)) = value(g(i));
        end
    end

    % Update iteration counter
    iter = iter + 1;
    if iter >= maxiter
        fprintf('Max iterations met.\n');
        break;
    end
    
    % Check termination conditions (relative improvement)
    rel_improve = (J_old - value(J))/abs(J_old);
    if rel_improve < epsilon/100
        fprintf('Relative improvement stalled.\n');
        break;
    else
        fprintf('J = %.5e, %.3f%% improvement.\n', value(J), 100*rel_improve);
        J_old = value(J);
    end
end
    
G = diag(G_k);

end

