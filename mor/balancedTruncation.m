function [A, B, C, W, V, S] = balancedTruncation(Af, Bf, Cf, k, continuous)
%[A, B, C, W, V, S] = balancedTruncation(Af, Bf, Cf, k, continuous)
%
%Perform balanced truncation on a stable continuous (discrete) LTI system, assuming
%the pair (Af, Bf) is controllable and (Af, Cf) is observable.
%Balanced truncation model order reduction to order k
%
% Inputs: 
%   Af, Bf, Cf: full order system matrices
%   k: reduced model order
%   continuous: true for continuous time system, false for discrete time
%
% Returns:
%   A, B, C: reduced order system matrices
%   W, V: Petrov-Galerkin projection matrices W'*V = I
%   S: Hankel singular values of balanced system

% Using the Balancing-free square-root algorithm, Antoulas Ch 7.3 pg 222
if continuous
    U = lyapchol(Af, Bf);
    L = lyapchol(Af', Cf');
else
    U = dlyapchol(Af, Bf);
    L = dlyapchol(Af', Cf');
end
U = U';
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

