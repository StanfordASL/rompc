function [EAdtBd] = recursiveMatMultiply(Ad, Bd, E, tau)
%[EAdtBd] = recursiveMatMultiply(Ad, Bd, E, tau)
%
% For the discrete time system x+ = Ad*x + Bd*u this computes the
% set of matrices [E*Bd, E*Ad*Bd, ... , E*Ad^tau-1*Bd] by recursively
% performing the matrix multiplications.
%
% Inputs:
%   Ad: discrete time system matrix
%   Bd: discrete time system matrix
%   E (optional): output matrix or [] to use I
%   tau: number of time steps for recursion
%
% Outputs:
%   EAdtBd: cell array [E*Bd, E*Ad*Bd, ... , E*Ad^tau-1*Bd]

fprintf('Using recursive matrix multiplications.\n');
if isempty(E)
    E = eye(size(Ad, 1));
end

fprintf('Precomputing some matrices.      ');
EAdtBd = cell(1, tau); % to get [E*Bd, E*Ad*Bd, ... , E*Ad^tau-1*Bd]
V = Bd;
EAdtBd{1} = full(E*Bd);
for i = 2:tau
    V = Ad*V;
    EAdtBd{i} = E*V;
    fprintf('\b\b\b\b\b%3.0f%% ', round(100*(i/(tau))));
end
fprintf('\n');



end

