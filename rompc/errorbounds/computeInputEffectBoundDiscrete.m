function [Delta_z, Delta_u] = computeInputEffectBoundDiscrete(ROM, ERROR, Z, U, Xbar, Wnoise, Vnoise, tau, opt)
%[Delta_z, Delta_u] = computeInputEffectBoundDiscrete(ROM, ERROR, Z, U, Xbar, Wnoise, Vnoise, tau, opt)
%
%Computes the bound Delta^(2) defined in Lorenzetti et al. 2020 for
%discrete time systems.
%
% Inputs:
%   ROM: struct containing discrete time ROM matrices (A, B, Bw, C, H)
%   ERROR: struct containing discrete time Ae, Be, Ge, Ez, and Eu error system matrices
%   Z: performance variable constraints (Polyhedron)
%   U: control constraints (Polyhedron)
%   Xbar: reduced state bounds (Polyhedron)
%   Wnoise: process noise bound (Polyhedron), or 0 if no noise
%   Vnoise: measurement noise bound (Polyhedron), or 0 if no noise
%   tau: backward time horizon
%   opt: optional arguments, including:
%       - solver: specifies the solver to use
%       - errormats_datapath: path to file where precomputed matrices are
%                             stored or loaded from
%
% Outputs:
%   Delta_z:
%   Delta_d:

n = size(ROM.A,1);
m = size(ROM.B,2);
p = size(ROM.C,1);
E = [ERROR.Ez; ERROR.Eu];

% Solver
if ~isfield(opt, 'solver')
    opt.solver = 'cplex';
end

% Check whether the matrices that are computed should be saved or loaded
if isfield(opt, 'errormats_datapath')

    % Load data matrices from file
    if exist(opt.errormats_datapath, 'file')
        tau_spec = tau;
        fprintf('******************************\n');
        fprintf('Loading error data matrices from %s.\nDelete this file for recomputation.\n', opt.errormats_datapath);
        fprintf('******************************\n\n');
        if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
            load(opt.errormats_datapath, 'EAetBe', 'EAetGe', 'tau');
        else
            load(opt.errormats_datapath, 'EAetBe', 'tau');
        end
        compute = false;

        % Check that tau matches with data
        if tau_spec <= tau
            tau = tau_spec;
        else
            fprintf('Using tau = %d because of stored error mat data.\n', tau);
        end

    % Compute data matrices and then save to file
    elseif ~exist(opt.errormats_datapath, 'file')
        fprintf('Will save data matrices to %s\n', opt.errormats_datapath);
        savedata = true;
        compute = true;
    end

% Compute data matrices and do not save
else
    savedata = false;
    compute = true;
end

% First compute useful matrices, this may take a while
if compute
    if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
        [EAetBe] = recursiveMatMultiply(ERROR.Ae, ERROR.Be, E, tau);
        [EAetGe] = recursiveMatMultiply(ERROR.Ae, ERROR.Ge, E, tau);
    else
        [EAetBe] = recursiveMatMultiply(ERROR.Ae, ERROR.Be, E, tau);
        EAetGe = [];
    end

    % Save the data
    if savedata
        fprintf('Saving data matrices to %s.\n\n', opt.errormats_datapath);
        save(opt.errormats_datapath, 'EAetBe', 'EAetGe', 'tau');
    end
end


% Define decision variables
ubar = {};
xbar = {};
fprintf('Creating decision variables.      ');
for i = 1:tau+tau % k = [-tau, ..., -1] and [0, ..., tau-1]
    ubar{i} = sdpvar(m, 1);
    xbar{i} = sdpvar(n, 1);
    fprintf('\b\b\b\b\b%3.0f%% ', round(100*(i/(tau+tau))));
end
fprintf('\n');

% Constraints
constraints = [];

fprintf('Building linear program.\n');
if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
    mw = size(Wnoise.A, 2);
    
    % Noise variables
    v = {};
    w = {};
    for i = 1:tau % k = [0, ..., tau-1]
        v{i} = sdpvar(p,1);
        w{i} = sdpvar(mw,1);
    end
    
    % Dynamics constraints
    for i = 1:tau+tau-1
        constraints = [constraints, xbar{i+1} - ROM.A*xbar{i} - ROM.B*ubar{i} == 0]; % dynamics of ROM
    end

    % Measurement constraints/process noise
    for i = 1:tau
        constraints = [constraints, Wnoise.A*w{i} <= Wnoise.b];
        constraints = [constraints, Vnoise.A*v{i} <= Vnoise.b];
    end
    
    % Cost function
    cost = 0;
    for i = 1:tau
        cost = cost + EAetBe{tau+1-i}*[xbar{tau+i}; ubar{tau+i}] + EAetGe{tau+1-i}*[w{i}; v{i}];
    end

% Otherwise do not include noise terms
else
    fprintf('No noise inputted, ignoring noise terms.\n');
    % Dynamics constraints
    for i = 1:tau+tau-1
        constraints = [constraints, xbar{i+1} - ROM.A*xbar{i} - ROM.B*ubar{i} == 0]; % dynamics of ROM
    end
    
    % Cost function
    cost = 0;
    for i = 1:tau
        cost = cost + EAetBe{tau+1-i}*[xbar{tau+i}; ubar{tau+i}];
    end
end

for i = 1:tau+tau
    constraints = [constraints, Xbar.A*xbar{i} <= Xbar.b]; % reduced state constraints
    constraints = [constraints, U.A*ubar{i} <= U.b]; % control constraints
    constraints = [constraints, Z.A*ROM.H*xbar{i} <= Z.b]; % reduced performance constraints
end


% Solve LPs to get p^Tx <= b
ops = sdpsettings('verbose', 0, 'solver', opt.solver, 'savesolveroutput', 1);
if strcmp(opt.solver, 'cplex')
    ops.cplex.lpmethod = 4;
    ops.cplex.barrier.convergetol = 1e-04;
end
n_z = size(Z.A,1);
Delta_z = zeros(n_z, 1);
fprintf('Computing bounds for Hf*e error.\n');
for i = 1:n_z
    fprintf(1, 'Computing %d of %d. ', i, n_z);
    J = -cost(i);
    diagnostics = optimize(constraints, J, ops);
    Delta_z(i) = -value(J);
    fprintf('Value = %.2f%%.\n', 100*Delta_z(i)/Z.b(i));

    % Cplex
    if diagnostics.problem ~= 0
        disp('Solver failed.');
        disp(diagnostics.info)
    end
end

n_u = size(U.A,1);
Delta_u = zeros(n_u, 1);
fprintf('Computing bounds for K*d error.\n');
for i = 1:n_u
    fprintf(1, 'Computing %d of %d. ', i, n_u);
    J = -cost(n_z + i);
    diagnostics = optimize(constraints, J, ops);
    Delta_u(i) = -value(J);
    fprintf('Value = %.2f%%.\n', 100*Delta_u(i)/U.b(i));

    % Cplex
    if diagnostics.problem ~= 0
        disp('Solver failed.');
        disp(diagnostics.info)
    end
end

end


