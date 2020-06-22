function [C] = computeC_VertexEnum(B, G, X, Y)
%[C] = computeC_VertexEnum(B, G, X, Y)
%
%Computes a tight bound on the norm ||B*[x, y]^T||_G <= C subject to the
%bounds on the variables x and y. Performs an exhaustive search
%over all vertices of bounding polyhedron.
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
fprintf('Computing vertex representations.\n');
X.computeVRep();
Y.computeVRep();
Nx = size(X.V, 1);
Ny = size(Y.V, 1);
N = Nx*Ny;
vals = zeros(N, 1);
idxs = round(1:(N/100):N);
fprintf('Checking %d vertices. ', N);
fprintf('  0%%');
for i = 1:Nx
    for j = 1:Ny
        idx = (i - 1)*Ny + j;
        r = [X.V(i, :), Y.V(j, :)];
        vals(idx) = sqrt(r*M*r');
        if ismember(idx, idxs)
            fprintf('\b\b\b\b%3.0f%%', round(100*(idx/N)));
        end
    end
end
C = max(vals);
fprintf('\nComputed C = %0.2f.\n', C);
end

