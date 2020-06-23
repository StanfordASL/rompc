function [C, T] = x2CT(x, params)
%[C, T] = x2CT(x, params)
%
%Converts the state vector x into the vectors for C [mol/L] and T [K] based
%on the normalization values in params as well as the linearization point

C_nondim = x(1:params.P.N);
T_nondim = x(params.P.N+1:end);
C = C_nondim*params.P.Cnorm + params.Cstar;
T = T_nondim*params.P.Tnorm + params.Tstar;
end

