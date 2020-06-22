function [A, B, C, W, V, S] = balancedTruncationDiscrete(Af, Bf, Cf, k)
%[A, B, C, W, V, S] = balancedTruncationDiscrete(Af, Bf, Cf, k)
%
% Perform balanced truncation on a stable discrete LTI system, assuming
% (Af, Bf) is controllable and (Af, Cf) is observable.
%   Balanced truncation model order reduction to order k
% Inputs: 
%       Af, Bf, Cf: full order system matrices
%       k: reduced model order
%
% Returns:
%       A, B, C: reduced order system matrices
%       W, V: Petrov-Galerkin projection matrices W'*V = I
%       S: singular values


% Using the Balancing-free square-root algorithm, Antoulas Ch 7.3 pg 222
U = dlyapchol(Af, Bf);
U = U';
L = dlyapchol(Af', Cf');
L = L';
[W, S, V] = svd(U'*L);
[X, ~] = qr(U*W);
[Y, ~] = qr(L*V);

Xk = X(:,1:k);
Yk = Y(:,1:k);
W = Yk*inv(Xk'*Yk);
V = Xk;

A = W'*Af*V;
B = W'*Bf;
C = Cf*V;
S = diag(S);

end

