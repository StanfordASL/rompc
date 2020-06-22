function [ERROR] = buildErrorSystem(FOM, ROM, CTRL, Z, U, Wnoise, Vnoise)
%[ERROR] = buildErrorSystem(FOM, ROM, CTRL, Z, U, Wnoise, Vnoise)
%
%Builds error system for continuous or discrete time FOM and ROM systems,
%along with a defined controller CTRL and optional polyhedral constraint
%sets Z and U.

nf = size(FOM.Af,1);
n = size(ROM.A,1);
m = size(FOM.Bf,2);
p = size(FOM.Cf,1);
o = size(FOM.Hf,1);
Pp = eye(nf) - ROM.V*inv(ROM.W'*ROM.V)*ROM.W';
Ae = [FOM.Af, FOM.Bf*CTRL.K; 
      CTRL.L*FOM.Cf, ROM.A + ROM.B*CTRL.K - CTRL.L*ROM.C];
Be = [Pp*FOM.Af*ROM.V, Pp*FOM.Bf; 
      sparse(n, n), sparse(n, m)];

if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
    mw = size(FOM.Bfw, 2);
    Ge = [FOM.Bfw, sparse(nf, p); sparse(n, mw), CTRL.L];
else
    Ge = [0];
end

if isa(Z, 'Polyhedron') && isa(U, 'Polyhedron')
    Ez = [Z.A*FOM.Hf*speye(nf), sparse(size(Z.A, 1), n)];
    Eu = [sparse(size(U.A, 1), nf), U.A*CTRL.K*speye(n)];
else
    Ez = 0;
    Eu = 0;
end

ERROR = struct('Ae', Ae, 'Be', Be, 'Ge', Ge, 'Ez', Ez, 'Eu', Eu);
end

