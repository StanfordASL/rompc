function [ROMPC, XF, SP, Zbar, Ubar, OPT] = buildROMPC(ROM, Z, U, EBOUND, N, opt)
%[ROMPC, XF, SP, Zbar, Ubar, OPT] = buildROMPC(ROM, Z, U, EBOUND, N, opt)
%
%Defines the ROMPC problem for a ROM with given constraints Z and U on the
%FOM and computed error bounds to tighten these constraints.
%
% Inputs:
%   ROM:
%   Z, U: polyhedral constraint sets
%   EBOUND: error structure
%   N: ROMPC time horizon (integer)
%   opt: optional arguments
%       - continuous (discrete): actually required
%       - eta_2tau: value for constraint tightening, default = 0
%       - setpoints: struct containing:
%               - FOM: full order model structure
%               - CTRL: controller data structure
%               - T: matrix specifying tracking variables r = Tz
%               - r: cell array of tracking column vectors, size(t,1)
%       - solver: solver for YALMIP (string), e.g. 'cplex' or 'mosek'
%       - output_x1: true have object output u(k) and x(k+1), false to just
%                    output u(k) (default)
%
% Returns:
%   ROMPC: without opt.setpoints returns ROMPC optimization object for the
%          origin. If opt.setpoints, then ROMPC is a cell array containing
%          optimization objects for origin + setpoints
%   Xf: without opt.setpoints returns terminal set for origin. Otherwise a
%       cell array like ROMPC
%   SP: without opt.setpoints is empty struct. Otherwise contains cell
%       array of structs that have the setpoints information
%   Zbar: tightened constraints
%   Ubar: tightened control constraints

if isfield(opt, 'continuous')
    continuous = opt.continuous;
elseif isfield(opt, 'discrete')
    continuous = ~opt.discrete;
else
    fprintf('Need to specify if the problem is discrete or continuous time.\n');
    return;
end

if ~isfield(opt, 'eta_2tau')
    fprintf('Ignoring eta_2tau term.\n');
    opt.eta_2tau = 0;
end

Nr = 0;
if isfield(opt, 'setpoints')
    if ~isfield(opt, 'FOM') || ~isfield(opt, 'CTRL')
        fprintf('Need to include opt.FOM and opt.CTRL to compute setpoints.\n');
    end
    
    % Extract matrix such that r = Tz
    if ~isfield(opt.setpoints, 'T')
        T = eye(o);
    else
        T = opt.setpoints.T;
    end
    
    if ~isfield(opt.setpoints, 'r') || (isfield(opt.setpoints, 'r') && size(r,1) ~= size(T,1))
        fprintf('Need valid setpoints in cell array opt.setpoints.r\n');
    else
        r = opt.setpoints.r;
    end
    Nr = length(r);
end

if continuous
    if ~isfield(EBOUND, 'dt')
        fprintf('EBOUND needs to specify time discretization.\n');
    end
    [Ad, Bd, ~] = zoh(EBOUND.dt, ROM.A, ROM.B, 0);
else
    Ad = ROM.A;
    Bd = ROM.B;
end

n = size(ROM.A,1);
m = size(ROM.B,2);
o = size(ROM.H,1);

% Tighten constraints based on error bounds
[~, ~, Zbar, Ubar] = tightenConstraints(Z, U, EBOUND, opt.eta_2tau);

% Compute a terminal controller and cost for MPC stability guarantees
[K, P, ~] = dlqr(Ad, -Bd, ROM.Q, ROM.R);

% Compute terminal set for origin
fprintf('Computing terminal set.\n');
Zf = Polyhedron(Z.A*ROM.H, Z.b);
Xbar1 = boundingBox(Zf, opt);
Xbar2 = Polyhedron(Zbar.A*ROM.H, Zbar.b);
Xbar = Xbar1.intersect(Xbar2);
Xf = computeTerminalSet(Ad, Bd, eye(n), K, Xbar, Ubar, zeros(n,1), zeros(m,1), opt);
Xf.minHRep();

% ROMPC for origin tracking
[ROMPC, OPT] = buildMPC(Ad, Bd, P, ROM.Q, ROM.R, ROM.H, Zbar, Ubar, Xf, N, zeros(n,1), zeros(m,1), opt);
SP = struct();
XF = Xf;

if Nr > 0
    ROMPC = {ROMPC};
    OPT = {OPT};
    XF = {Xf};
    SP = {SP};
    for i = 1:Nr
        fprintf('Computing setpoint.\n');
        if continuous
            [sp.xfss, sp.uss, sp.xbarss, sp.ubarss, sp.xhatss] = computeSteadyStateContinuousTime(opt.FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r{i});
        else
            [sp.xfss, sp.uss, sp.xbarss, sp.ubarss, sp.xhatss] = computeSteadyStateDiscreteTime(opt.FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r{i});
        end
        fprintf('Computing terminal set.\n');
        Xf = computeTerminalSet(Ad, Bd, eye(n), K, Xbar, Ubar, sp.xbarss, sp.ubarss, opt);
        Xf.minHRep();
        
        % Set up ROMPC object
        [ROMPC{i+1}, OPT{i+1}] = buildMPC(Ad, Bd, P, ROM.Q, ROM.R, ROM.H, Zbar, Ubar, Xf, N, sp.xbarss, sp.ubarss, opt);
        
        % Store new data
        XF{i+1} = Xf;
        SP{i+1} = sp;
    end
end

end

