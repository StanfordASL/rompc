function X = boundingBox(Z, opt)
%X = boundingBox(Z, opt)
%
%Computes a bounding box around the set Z. If Z is not bounded the bounding
%box may have some terms be infinity.
%
% Inputs:
%   Z: Polyhedron object defining the set
%   opt: opt.solver specifies which solver to use
%
% Returns:
%   X: bounding box Polyhedron

N = Z.Dim;
M = size(Z.b, 1);

[P] = rectanglePolytope(zeros(1,N), zeros(1,N));

% Use the rectangle vectors and solve the support functions with them
[prob] = supportFunction(N, M, opt);
hZ = prob([], Z.A, Z.b);
b_x = zeros(2*N,1);
for i = 1:2*N
    [~,~,~,~,~,solver] = hZ(P.A(i, :)');
    if solver.problem == 2
        b_x(i) = inf;
    else
        b_x(i) = -solver.solveroutput.fval;
    end
end
X = Polyhedron(P.A, b_x);
end

