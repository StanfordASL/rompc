function [ROMPC, Xf, SP, Zbar, Ubar] = buildROMPC(ROM, Z, U, EBOUND, N, opt)
%[ROMPC, Xf, SP, Zbar, Ubar] = buildROMPC(ROM, Z, U, EBOUND, N, opt)
%
%Defines the ROMPC problem for a ROM with given constraints Z and U on the
%FOM and computed error bounds to tighten these constraints.
%
% Inputs:
%   ROM: reduced order model structure
%   Z, U: polyhedral constraint sets
%   EBOUND: error structure
%   N: ROMPC time horizon (integer)
%   opt: options arguments
%       - continuous (discrete): true or false
%       - eta_2tau: value for constraint tightening, default = 0
%       - setpoint: struct containing:
%               - FOM: full order model structure
%               - CTRL: controller data structure
%               - T: matrix specifying tracking variables r = Tz
%               - r: cell array of tracking column vectors, size(t,1)
%                   defaults to the origin if opt.setpoint is not defined
%       - solver: solver for YALMIP (string), e.g. 'cplex' or 'mosek'
%       - xf_datapath: path for where to save or load terminal set
%
% Returns:
%   ROMPC: struct with 
%       - rompc: precompiled YALMIP problem
%       - prob: uncompiled YALMIP optimization problem
%       - dt: time discretization for ROMPC, -1 for discrete time problem
%   Xf: terminal set for defined opt.setpoint, or the origin
%   SP: setpoint steady state information
%   Zbar: tightened constraints
%   Ubar: tightened control constraints

n = size(ROM.A,1);
m = size(ROM.B,2);
o = size(ROM.H,1);

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

% Tighten constraints based on error bounds
[~, ~, Zbar, Ubar] = tightenConstraints(Z, U, EBOUND, opt.eta_2tau);

% Handle setpoints
if isfield(opt, 'setpoint')
    if ~isfield(opt.setpoint, 'FOM') || ~isfield(opt.setpoint, 'CTRL')
        error('Need to include opt.setpoint.FOM and opt.setpoint.CTRL to compute setpoints.');
    end
    FOM = opt.setpoint.FOM;
    CTRL = opt.setpoint.CTRL;
    
    % Remove them because no need to have them saved later
    opt.setpoint = rmfield(opt.setpoint, 'FOM');
    opt.setpoint = rmfield(opt.setpoint, 'CTRL');
    
    % Extract matrix such that r = Tz
    if ~isfield(opt.setpoint, 'T')
        T = eye(o);
    else
        T = opt.setpoint.T;
    end
    
    if ~isfield(opt.setpoint, 'r') || (isfield(opt.setpoint, 'r') && ...
                size(opt.setpoint.r,1) ~= size(opt.setpoint.T,1))
        error('Need valid setpoint in opt.setpoint.r');
    else
        r = opt.setpoint.r;
    end
    fprintf('Computing setpoint.\n');
    if continuous
        [SP.xfss, SP.uss, SP.xbarss, SP.ubarss, SP.xhatss] = computeSteadyStateContinuousTime(FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r);
    else
        [SP.xfss, SP.uss, SP.xbarss, SP.ubarss, SP.xhatss] = computeSteadyStateDiscreteTime(FOM, ROM, CTRL, Z, U, Zbar, Ubar, T, r);
    end
else
    SP.xbarss = zeros(n,1);
    SP.ubarss = zeros(m,1);
end

% Discretize ROM if continuous time
if continuous
    if ~isfield(EBOUND, 'dt')
        fprintf('EBOUND needs to specify time discretization.\n');
    end
    [Ad, Bd, ~] = zoh(EBOUND.dt, ROM.A, ROM.B, 0);
else
    Ad = ROM.A;
    Bd = ROM.B;
end

% Compute terminal set
if isfield(opt, 'xf_datapath') && exist(opt.xf_datapath, 'file')
    fprintf('******************************\n');
    fprintf('Loading terminal set from %s.\nDelete this file to recompute.\n', opt.xf_datapath);
    fprintf('******************************\n\n');
    load(opt.xf_datapath, 'Xf', 'K', 'P', 'SP');
    fprintf('Using setpoint that was used to define saved terminal set.\n');
else
    % Compute a terminal controller and cost for MPC stability guarantees
    [K, P, ~] = dlqr(Ad, -Bd, ROM.Q, ROM.R);

    % Compute terminal set for origin
    fprintf('Computing terminal set.\n');
    Zf = Polyhedron(Z.A*ROM.H, Z.b);
    Xbar1 = boundingBox(Zf, opt);
    Xbar2 = Polyhedron(Zbar.A*ROM.H, Zbar.b);
    Xbar = Xbar1.intersect(Xbar2);
    Xf = computeTerminalSet(Ad, Bd, eye(n), K, Xbar, Ubar, SP.xbarss, SP.ubarss, opt);
    Xf.minHRep();
    
    if isfield(opt, 'xf_datapath')
        fprintf('Saving terminal set data to %s\n', opt.xf_datapath);
        save(opt.xf_datapath, 'Xf','K','P','SP');
    end
end

% ROMPC for setpoint tracking
[rompc, prob] = buildMPC(Ad, Bd, P, ROM.Q, ROM.R, ROM.H, Zbar, Ubar, Xf, N, SP.xbarss, SP.ubarss, opt);
ROMPC.rompc = rompc;
ROMPC.prob = prob;
if continuous
    ROMPC.dt = EBOUND.dt;
else
    ROMPC.dt = -1;
end

end
