function [norm] = matrix2Norm(A)
%[norm] = matrix2Norm(A)
%
%Computes the matrix 2 norm of A, ||A||_2 via eigenvalue decomposition.

[~, lam] = eig(A'*A);
norm = max(sqrt(diag(lam)));

end

