function [xfss, uss, xbarss, ubarss, xhatss] = computeSteadyStateDiscreteTime(FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r)
%[xfss, uss, xbarss, ubarss, xhatss] = computeSteadyStateDiscreteTime(FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r)
%
%   Computes the steady state values for the full order system and reduced
%   order target states based on estimator and controller for a discrete
%   time system.
%
% Inputs:
%   FOM: struct containing continuous time FOM matrices (Af, Bf, Bfw, Cf, Hf)
%   ROM: struct containing continuous time ROM matrices (A, B, Bw, C, H, V, W)
%   CTRL: controllers K and L
%   Z: performance variable constraints (Polyhedron)
%   U: control constraints (Polyhedron)
%   Zbar: tightened performance variable constraints (Polyhedron)
%   Ubar: tightened control constraints (Polyhedron)
%   T: z_track = T*z
%   r: z_track = r

nf = size(FOM.Af,1); % full state size
m = size(FOM.Bf,2); % control size
n = size(ROM.A, 1); % reduced order state size
t = size(T, 1); % tracking variable dimension

% Compute the steady state values
Ass = [FOM.Af - eye(nf), FOM.Bf;
       T*FOM.Hf, zeros(t, m)];
bss = [zeros(nf,1);
       r];
sol = pinv(full(Ass))*bss;
xfss = sol(1:nf);
uss = sol(nf+1:end);

% Compute reduced order steady state value
xhatss = linsolve(eye(n) - (ROM.A - CTRL.L*ROM.C), ROM.B*uss + CTRL.L*FOM.Cf*xfss);
sol = linsolve([eye(n) - ROM.A, -ROM.B; -CTRL.K, eye(m)],[zeros(n, 1); uss - CTRL.K*xhatss]);
xbarss = sol(1:n);
ubarss = sol(n+1:end);

% Check if these values are within the tightened bounds
if Z.contains(FOM.Hf*xfss) && U.contains(uss)
    fprintf('Full-order steady state state and control are valid.\n');
else
    fprintf('The tracking point is not feasible with respect to original constraints, returning the origin.\n');
    xfss = zeros(nf,1);
    uss = zeros(m,1);
    xbarss = zeros(n,1);
    ubarss = uss;
    xhatss = xbarss;
end

% Check if these values are within the tightened bounds
if Zbar.contains(ROM.H*xbarss) && Ubar.contains(ubarss)
    fprintf('Reduced-order steady state state and control are valid\n');
else
    fprintf('The tracking point is not feasible with respect to reduced (tightened) constraints, returning the origin.\n');
    xfss = zeros(nf,1);
    uss = zeros(m,1);
    xbarss = zeros(n,1);
    ubarss = uss;
    xhatss = xbarss;
end

end

