function [C] = computeC_HyperRect(B, G, X, Y)
%[C] = computeC_HyperRect(B, G, X, Y)
%
%Computes a tight bound on the norm ||B*[x, y]^T||_G <= C subject to the
%bounds on the variables x and y. Performs an exhaustive search
%over all vertices of bounding polyhedron, assuming they are
%hyperrectangles.
%
% Inputs:
%   B: matrix
%   G: positive definite weighting matrix
%   X: bounding set for x (compact)
%   Y: bounding set for y (compact)
%
% Outputs:
%   C: norm bound

M = B'*G*B;

% Check each vertex
fprintf('Checking vertices. ');
n = size(X.A, 2);
m = size(Y.A, 2);

fprintf('  0%%');
[C, N] = checkVertices(0, zeros(n,1), zeros(m,1), M, X, Y, n, m, 1, 1, 0);
fprintf('\nEvaluated %d vertices.\n', N);
fprintf('Computed C = %0.2f.\n', C);
end


function [val_max, N] = checkVertices(val_max, x, u, M, X, Y, n, m, ix, iu, N)
% If x is not built keep going
if ix <= n
    % Look for positive coordinate
    x_cur = x;
    x_cur(ix) = X.b(X.A(:, ix) == 1);
    [val_max, N] = checkVertices(val_max, x_cur, u, M, X, Y, n, m, ix + 1, iu, N);

    % Look for negative coordinate
    x_cur = x;
    x_cur(ix) = -X.b(X.A(:, ix) == -1);
    [val_max, N] = checkVertices(val_max, x_cur, u, M, X, Y, n, m, ix + 1, iu, N);
    
else % the full vector x has been built
    if iu <= m
        % Look for positive coordinate
        u_cur = u;
        u_cur(iu) = Y.b(Y.A(:, iu) == 1);
        [val_max, N] = checkVertices(val_max, x, u_cur, M, X, Y, n, m, ix, iu + 1, N);

        % Look for negative coordinate
        u_cur = u;
        u_cur(iu) = -Y.b(Y.A(:, iu) == -1);
        [val_max, N] = checkVertices(val_max, x, u_cur, M, X, Y, n, m, ix, iu + 1, N);
    
    else % the full vector has been built
        r = [x; u];
        val = sqrt(r'*M*r);
        if val > val_max
            val_max = val;
        end
        N = N + 1;
        
        if mod(100*N/(2^(n+m)), 1) == 0
           fprintf('\b\b\b\b%3.0f%%', round(100*(N/2^(n+m))));
        end
    end
end

end

