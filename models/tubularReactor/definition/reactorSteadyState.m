function [Cstar, Tstar] = reactorSteadyState(TJ1, TJ2, TJ3, P)
%[Cstar, Tstar] = reactorSteadyState(TJ1, TJ2, TJ3, P)
%
%Computes a steady state for the reactor given the steady state jacket
%temperatures and the model parameters.

% Compute steady state operating point
    
Cstar = zeros(1,P.N);
Tstar = zeros(1,P.N);
Cstar(1) = P.Cin/P.Cnorm;
Tstar(1) = P.Tin/P.Tnorm;
for i = 2:P.N
    if i <= P.N/3
        u = TJ1/P.Tnorm;
    elseif i > P.N/3 && i <= 2*P.N/3
        u = TJ2/P.Tnorm;
    else
        u = TJ3/P.Tnorm;
    end
    Cstar(i) = Cstar(i-1)*(1 - P.k0*P.dz*exp(-P.E/(P.R*P.Tnorm*Tstar(i-1)))/P.v);
    Tstar(i) = Tstar(i-1)*(1 - P.Hr*P.dz/P.v) ...
               + (P.Gr*P.dz*P.Cnorm/(P.v*P.Tnorm))*Cstar(i-1)*exp(-P.E/(P.R*P.Tnorm*Tstar(i-1))) ...
               + (P.Hr*P.dz/P.v)*u;
end
Tstar = Tstar'*P.Tnorm;
Cstar = Cstar'*P.Cnorm';
end
