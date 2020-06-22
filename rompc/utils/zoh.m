function [Ad, Bd, Bwd] = zoh(dt, A, B, Bw)
%[Ad, Bd, Bwd] = zoh(dt, A, B, Bw)
%
%Discretizes a continuous time system by assuming a zero order hold on the
%system inputs. Can take two input matrices, B and Bw, for control and
%disturbance respectively. Computes discretization by solving for matrix
%exponential.
%
% Inputs:
%   dt: discretization time
%   A: system matrix
%   B: control matrix
%   Bw: disturbance matrix, or 0 if not included
%   
% Returns:
%   Ad: ZOH time system matrix
%   Bd: ZOH control matrix
%   Bwd: ZOH disturbance matrix

n = size(A,1);
m = size(B,2);
mw = size(Bw,2);

if size(Bw,1) == n
    ZOH = expm(dt*[A, B, Bw; zeros(m+mw, n), zeros(m+mw)]);
    Ad = ZOH(1:n,1:n);
    Bd = ZOH(1:n, n+1:n+m);
    Bwd = ZOH(1:n, n+m+1:end);
else
    ZOH = expm(dt*[A, B; zeros(m,n), zeros(m)]);
    Ad = ZOH(1:n,1:n);
    Bd = ZOH(1:n, n+1:end);
    Bwd = 0;
end

end

