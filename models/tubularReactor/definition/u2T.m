function [uT] = u2T(u, params)
%[uT] = u2T(u, params)
%
%Converts the control vector u into the thermal jacket temperatures based
%on the normalization values in params as well as the linearization point

uT = u*params.P.Tnorm + params.ustar;
end

