function [P] = rectanglePolytope(ub, lb)
%Creates a MPT3 Polyhedron object that is a rectangle centered around
%origin that has specified upper and lower bounds
% Inputs:
%   ub: array of upper bounds for each dimension
%   lb: array of lower bounds for each dimension
%
% Outputs:
%   P: rectangular polyhedron

n = length(ub);
A = [];
b = [];
basis = eye(n);
for i = 1:n
    A = [A; basis(i,:); -basis(i,:)];
    b = [b; ub(i); -lb(i)];
end

P = Polyhedron(A, b);
end

