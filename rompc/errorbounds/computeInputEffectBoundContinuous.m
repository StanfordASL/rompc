function [Delta_z, Delta_u] = computeInputEffectBoundContinuous(ROM, ERROR, Z, U, Xbar, Wnoise, Vnoise, tau, opt)
%[Delta_z, Delta_u] = computeInputEffectBoundContinuous(ROM, ERROR, Z, U, Xbar, Wnoise, Vnoise, tau, opt)
%
%Computes the bound Delta^(2) defined in Lorenzetti et al. 2020 but for
%continuous time systems instead of discrete time systems.
%
% Inputs:
%   ROM: struct containing continuous time ROM matrices (A, B, Bw, C, H)
%   ERROR: struct containing continuous time Ae, Be, Ge, Ez, and Eu error system matrices
%   Z: performance variable constraints (Polyhedron)
%   U: control constraints (Polyhedron)
%   Xbar: reduced state bounds (Polyhedron)
%   Wnoise: process noise bound (Polyhedron), or 0 if no noise
%   Vnoise: measurement noise bound (Polyhedron), or 0 if no noise
%   tau: backward time horizon
%   opt: optional arguments, including:
%       - errormats_datapath: path to file where precomputed matrices are
%                             stored to or loaded from
%       - dt: time discretization
%       - int_method: 'trapezoid' or 'leftRiemann'
%       - mat_exp_method: 'expm', 'expmv', or 'implicitEuler'
%
% Outputs:
%   Delta_z:
%   Delta_d:

n = size(ROM.A,1);
m = size(ROM.B,2);
p = size(ROM.C,1);

% Integration scheme
if isfield(opt, 'int_method')
    int_method = opt.int_method;
else
    int_method = 'trapezoid';
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
            load(opt.errormats_datapath, 'EexpAetBe', 'EexpAetGe', 'dt', 'tau', 'N');
        else
            load(opt.errormats_datapath, 'EexpAetBe', 'dt', 'tau', 'N');
        end
        compute = false;

        % Check that dt matches data dt
        fprintf('Using data value of dt = %0.3f.\n', dt);
        
        % Check that tau matches with data
        if tau_spec <= tau
            N = floor(tau_spec/dt);
            tau = N*dt;
        else
            fprintf('Using tau = %d because of stored error mat data.\n', tau);
        end

        % Check that tau matches with data
        fprintf('Using data values of tau = %d and N = %d.\n', tau, N);       

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
    E = [ERROR.Ez; ERROR.Eu];

    if isfield(opt, 'mat_exp_method')
        method = opt.mat_exp_method;
    else
        fprintf('Computing data matrices by first computing matrix exponential.\n');
        method = 'expm';
    end

    if isfield(opt, 'dt')
        dt = opt.dt;
    else
        fprintf('opt.dt not defined, using dt = 0.01\n');
        dt = 0.01;
    end
    N = floor(tau/dt);
    tau = N*dt;
    fprintf('Using tau = %.3f and N = %d.\n', tau, N);

    if strcmp(method, 'expm')
        [Aed, ~, ~] = zoh(dt, ERROR.Ae, ERROR.Be, 0);
        if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
            [EexpAetBe] = recursiveMatMultiply(Aed, ERROR.Be, E, N+1);
            [EexpAetGe] = recursiveMatMultiply(Aed, ERROR.Ge, E, N+1);
        else
            [EexpAetBe] = recursiveMatMultiply(Aed, ERROR.Be, E, N+1);
            EexpAetGe = [];
        end
        
    elseif strcmp(method, 'expmv')
        if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
            [EexpAetBe] = recursiveMatExpAction(ERROR.Ae, ERROR.Be, E, dt, N);
            [EexpAetGe] = recursiveMatExpAction(ERROR.Ae, ERROR.Ge, E, dt, N);
        else
            [EexpAetBe] = recursiveMatExpAction(ERROR.Ae, ERROR.Be, E, dt, N);
            EexpAetGe = [];
        end
        
    elseif strcmp(method, 'implicitEuler')
        if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
            [EexpAetBe] = recursiveImplicitEuler(ERROR.Ae, ERROR.Be, E, dt, N);
            [EexpAetGe] = recursiveImplicitEuler(ERROR.Ae, ERROR.Ge, E, dt, N);
        else
            [EexpAetBe] = recursiveImplicitEuler(ERROR.Ae, ERROR.Be, E, dt, N);
            EexpAetGe = [];
        end
        
    else
        fprintf('Method %s not known.\n', method);
        return;
    end

    % Save the data
    if savedata
        fprintf('Saving data matrices to %s.\n\n', opt.errormats_datapath);
        save(opt.errormats_datapath, 'EexpAetBe', 'EexpAetGe', 'dt', 'tau', 'N');
    end
end


% Define decision variables
ubar = {};
xbar = {};
fprintf('Creating decision variables.      ');
for i = 1:N+N+1 % k = [-tau, ..., -dt] and [0, ..., tau]
    ubar{i} = sdpvar(m, 1);
    xbar{i} = sdpvar(n, 1);
    fprintf('\b\b\b\b\b%3.0f%% ', round(100*(i/(N+N+1))));
end
fprintf('\n');

% Discretize ROM
[Ad, Bd, ~] = zoh(dt, ROM.A, ROM.B, 0);

% Define quadrature weights
weights = NewtonCotes(int_method, N, dt);

% Constraints
constraints = [];

fprintf('Building linear program.\n');
if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
    mw = size(Wnoise.A, 2);
    
    % Noise variables
    v = {};
    w = {};
    for i = 1:N+1 % k = [0, dt, ..., tau]
        v{i} = sdpvar(p,1);
        w{i} = sdpvar(mw,1);
    end
    
    % Dynamics constraints
    for i = 1:N+N
        constraints = [constraints, xbar{i+1} - Ad*xbar{i} - Bd*ubar{i} == 0]; % dynamics of ROM
    end

    % Measurement constraints/process noise
    for i = 1:N+1
        constraints = [constraints, Wnoise.A*w{i} <= Wnoise.b];
        constraints = [constraints, Vnoise.A*v{i} <= Vnoise.b];
    end
    
    % Cost function
    cost = 0;
    for i = 1:N+1
        cost = cost + weights(i).*(EexpAetBe{N+2-i}*[xbar{N+i}; ubar{N+i}] + EexpAetGe{N+2-i}*[w{i}; v{i}]);
    end

% Otherwise do not include noise terms
else
    fprintf('No noise inputted, ignoring noise terms.\n');
    % Dynamics constraints
    for i = 1:N+N
        constraints = [constraints, xbar{i+1} - Ad*xbar{i} - Bd*ubar{i} == 0]; % dynamics of ROM
    end
    
    % Cost function
    cost = 0;
    for i = 1:N+1
        cost = cost + weights(i).*(EexpAetBe{N+2-i}*[xbar{N+i}; ubar{N+i}]);
    end
end

for i = 1:N+N+1
    constraints = [constraints, Xbar.A*xbar{i} <= Xbar.b]; % reduced state constraints
    constraints = [constraints, U.A*ubar{i} <= U.b]; % control constraints
    constraints = [constraints, Z.A*ROM.H*xbar{i} <= Z.b]; % reduced performance constraints
end

% Solve LPs to get p^Tx <= b
ops = sdpsettings('verbose', 0, 'solver', 'cplex', 'savesolveroutput', 1);
ops.cplex.lpmethod = 4;
ops.cplex.barrier.convergetol = 1e-04;
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
        Delta_z(i) = inf;
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
        Delta_u(i) = inf;
    end
end

end


