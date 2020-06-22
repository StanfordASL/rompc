function [Cbound] = computeC_UpperBound(B, G, X, Y)
%[Cbound] = computeC_UpperBound(B, G, X, Y)
%
%Computes a conservative upper bound on the norm ||B*[x, y]^T||_G
%subject to the bounds on the variables x and y. This value is
%an upper bound on the value C.
% Inputs:
%   B: matrix
%   G: positive definite weighting matrix
%   X: bounding set for x
%   Y: bounding set for y
%
% Outputs:
%   Cbound: norm bound

M = B'*G*B;
n = size(X.A, 2);
m = size(Y.A, 2);
A = [X.A, zeros(size(X.A,1),m); zeros(size(Y.A,1),n), Y.A];
b = [X.b; Y.b];

% Start by solving LP for gamma
x = sdpvar(size(A,2),1);
y = sdpvar(size(A,1),1); 
constraints = [A*x + y == b, y >= 0];
J = -sum(y);
ops = sdpsettings('verbose', 0, 'solver', 'cplex', 'savesolveroutput', 1);
diagnostics = optimize(constraints, J, ops);
gamma = -value(J);
if diagnostics.problem ~= 0
    disp('Solver failed.\n');
    disp(diagnostics.info)
    gamma = inf;
end

% Bound C
B = inv(A'*A)*A';
d = B*b;
vals = [d'*M*d];
for i = 1:size(A,1)
    ei = zeros(size(A,1),1); ei(i) = 1;
    v = d - gamma*B*ei;
    vals = [vals, v'*M*v];
end
Cbound = max(vals);
fprintf('Computed C upper bound = %0.2f.\n', Cbound);
end

