function [w] = NewtonCotes(method, N, dt)
%[w] = NewtonCotes(method, N, dt)
%
%Specifies the quadrature weights for different Newton-Coates
%schemes.
%
% Inputs: 
%   - method: either 'trapezoid' or 'leftRiemann'
%   - N: number of intervals to use
%   - dt: width of each interval
%
% Returns:
%   - w: quadrature weights

if strcmp(method, 'trapezoid')
    fprintf('Using trapezoid quadrature scheme.\n');
    w = dt*ones(1, N+1);
    w(1) = dt/2;
    w(N+1) = dt/2;
    
elseif strcmp(method, 'leftRiemann')
    fprintf('Using left Riemann quadrature scheme.\n');
    w = dt*ones(1, N+1);
    w(N+1) = 0;
else
    fprintf('Method %s not implemented.\n', method);
end
    
end


