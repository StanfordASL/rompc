function [Xf] = computeTerminalSet(A, B, H, K, Zbar, Ubar, xbarss, ubarss, opt)
%[Xf] = computeTerminalSet(A, B, H, Zbar, Ubar, xbarss, ubarss)
%
%Computes the terminal set Xf centered around xbarss. A and B should be for
%a discrete time system.
%
% Inputs:
%   A: reduced order system matrix
%   B: reduced order control matrix
%   H: reduced order performance variable matrix
%   K: controller matrix A+BK stable
%   Zbar: performance variable constraints (Polyhedron)
%   Ubar: control constraints (Polyhedron)
%   xbarss: steady state xbar
%   ubarss: steady state ubar
%   opt: opt.solver specifies solver to use

% Shift the sets to the origin
DZbar = Zbar - H*xbarss;
DUbar = Ubar - ubarss;

% Start out with the largest set that satisfies both state and control
% constraints
Delta = Polyhedron([DZbar.A*H; DUbar.A*K], [DZbar.b; DUbar.b]);
Delta.minHRep();

if ~Delta.isBounded()
    fprintf('Largest possible terminal set is not bounded.\n');
end

% Compute invariant set
AK = A + B*K;
fprintf('Computing invariant terminal set. \n');
[D, ~, ~] = maximalPI(AK, eye(size(AK,1)), Delta, opt);

% Convert the invariant set about the origin to actual Xf
Xf = D + xbarss;

end

