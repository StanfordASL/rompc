function [u] = T2u(uT, params)
%[u] = T2u(uT, params)
%
%Converts the control vector uT (thermal jacket temperatures) to normalized u based
%on the normalization values in params as well as the linearization point

u = (uT - params.ustar)/params.P.Tnorm;
end

