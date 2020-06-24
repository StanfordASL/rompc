function [A, B, C, W, V, S] = balancedTruncationContinuousUnstable(Af, Bf, Cf, k)
%[A, B, C, W, V, S] = balancedTruncationContinuousUnstable(Af, Bf, Cf, k)
%
% Perform balanced truncation on a unstable continuous LTI system, assuming
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

nf = size(Af,1);
m = size(Bf,2);
sysf = ss(full(Af), full(Bf), eye(nf), zeros(nf,m));
fprintf('Performing stable unstable decomposition.\n');
[S, NS] = stabsep(sysf); % decompose system into stable and unstable parts

% Use decomposed system to form stable subsystem
Afstable = S.A; Bfstable = S.B; Tfstable = S.C;
nfstable = size(Afstable,1);

% Use decomposed system to form unstable subsystem
Tfunstable = NS.C;
nfunstable = nf - nfstable;

T = [Tfunstable, Tfstable]; % such that T^-1 Af T = diag(Afunstable, Afstable)

% Compute reduced order model of the stable part of the system
nROM = k - nfunstable;
fprintf('Performing model reduction.\n');
Cftransformed = Cf*T;
Cfstable = Cftransformed(:, nfunstable + 1:end);
[~, ~, ~, Wstable, Vstable, S] = balancedTruncationContinuous(Afstable, Bfstable, Cfstable, nROM);

% Combine to make full ROM
W = ([eye(nfunstable), zeros(nfunstable, nROM);
          zeros(nfstable, nfunstable), Wstable]'*inv(T))';
V = T*[eye(nfunstable), zeros(nfunstable, nROM);
       zeros(nfstable, nfunstable), Vstable]; % A = W'*Af*V
A = W'*Af*V;
B = W'*Bf;
C = Cf*V;

end

