function [x] = CT2x(C, T, params)
%[x] = CT2x(C, T, params)
%
%Converts vectors for C [mol/L] and T [K] into the state vector x based
%on the normalization values in params as well as the linearization point

C_nondim = (C - params.Cstar)/params.P.Cnorm;
T_nondim = (T - params.Tstar)/params.P.Tnorm;
x = [C_nondim; T_nondim];
end

