function [Dz, Du, Zbar, Ubar] = tightenConstraints(Z, U, EBOUND, eta_2tau)
%[Dz, Du, Zbar, Ubar] = tightenConstraints(Z, U, EBOUND, eta_2tau)
%
% Computes the tightened constraints given the precomputed error bound
% parameters, and the bound eta_2tau. This is as defined in Lorenzetti et
% al. 2020. 

if isfield(EBOUND, 'eta_z') && isfield(EBOUND, 'eta_u') && isfield(EBOUND, 'D_1a') && isfield(EBOUND, 'D_1b')
    Dz = EBOUND.eta_z*(EBOUND.D_1a*eta_2tau + EBOUND.D_1b) + EBOUND.Dz_2;
    Du = EBOUND.eta_u*(EBOUND.D_1a*eta_2tau + EBOUND.D_1b) + EBOUND.Du_2;
else
    fprintf('Using only Dz_2 terms from error bound.\n');
    Dz = EBOUND.Dz_2;
    Du = EBOUND.Du_2;
end

% Tighten constraints
Zbar = Polyhedron(Z.A, Z.b - Dz);
Ubar = Polyhedron(U.A, U.b - Du);
end

