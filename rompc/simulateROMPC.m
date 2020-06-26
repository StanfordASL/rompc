function [DATA] = simulateROMPC(FOM, ROM, CTRL, ROM_CTRL, Zbar, T, X0, opt)
%[DATA] = simulateROMPC(FOM, ROM, CTRL, ROMPC, Zbar, T, X0, opt)
%
%Simulate the FOM with ROMPC control.
%
% Inputs:
%   FOM: struct containing FOM matrices (Af, Bf, Bfw, Cf, Hf)
%   ROM: struct containing ROM matrices (A, B, Bw, C, H, V, W)
%   CTRL: struct containing gains K and L
%   ROM_CTRL: controller for ROM, either is the precompiled ROMPC problem, a
%          structure containing the ROMPC problem information, or a gain
%          matrix for state feedback control
%   Zbar: tightened constraints used in ROMPC
%   T: time horizon, time steps for discrete and time (seconds, minutes)
%      for continuous time problems
%   X0: struct containing xf, xhat, xbar initial conditions
%   opt: options arguments, including
%       - continuous (discrete): true or false
%       - NOISE: struct containing W or V Polyhedrons for random bounded
%                disturbances, or wf/v matrices of size mw/p by T to
%                directly specify the noise terms at each time step
%       - sim_dt: for continuous time simulation, time step for integration
%       - rompc_dt: for continuous time simulation, discretization time
%                   used in ROMPC definition (ROMPC.dt)
%       - save_xf: true to save full state, default = false
%
% Returns:
%   DATA: struct containing simulation data as well as options and
%         controller query times

nf = size(FOM.Af, 1);
m = size(FOM.Bf, 2);
p = size(FOM.Cf, 1);
o = size(FOM.Hf, 1);

if ~isfield(opt, 'save_xf')
    opt.save_xf = false;
end

% Continuous or discrete time
if isfield(opt, 'continuous')
    continuous = opt.continuous;
    if ~isfield(opt, 'sim_dt') || ~isfield(opt, 'rompc_dt')
        error('Need to specify opt.sim_dt and opt.rompc_dt for continuous time problems');
    end
    
    % Convert time to number of time steps
    t = 0:opt.sim_dt:T;
    T = length(t);
    
    % How many steps before recomputing ROMPC
    update = opt.rompc_dt/opt.sim_dt;
    if floor(update) ~= update
        error('opt.rompc_dt should be integer multiple of opt.sim_dt');
    end
    
    % Discretize ROM and estimator for easy simulation
    [Ad, Bd, ~] = zoh(opt.sim_dt, ROM.A, ROM.B, 0);
    [Ade, Bde, Bwe] = zoh(opt.sim_dt, ROM.A - CTRL.L*ROM.C, ROM.B, CTRL.L);
    
    % Define preconditioner for backward Euler simulation
    S = speye(nf) - opt.sim_dt*FOM.Af;
    setup.type = 'nofill';
    [M1, M2] = ilu(sparse(S), setup);
    
elseif isfield(opt, 'discrete')
    continuous = ~opt.discrete;
    update = -1;
else
    error('Need to specify if the problem is discrete or continuous time.');
end

% Check what kind of controller was used for ROM
if isstruct(ROM_CTRL)
    precompiled = false;
    rompc = true;
elseif isa(ROM_CTRL, 'optimizer')
    precompiled = true;
    rompc = true;
elseif ismatrix(ROM_CTRL)
    rompc = false;
end

if ~isfield(opt, 'NOISE')
    opt.NOISE = struct();
    FOM.Bfw = zeros(nf,1);
else
    if isfield(opt.NOISE, 'W')
        mw = size(opt.NOISE.W.A,2);
        idxUB = find(opt.NOISE.W.A*ones(mw,1) == 1);
        idxLB = find(opt.NOISE.W.A*ones(mw,1) == -1);
        wfUB = opt.NOISE.W.b(idxUB);
        wfLB = -opt.NOISE.W.b(idxLB);
    end
    if isfield(opt.NOISE, 'V')
        idxUB = find(opt.NOISE.V.A*ones(p,1) == 1);
        idxLB = find(opt.NOISE.V.A*ones(p,1) == -1);
        vUB = opt.NOISE.V.b(idxUB);
        vLB = -opt.NOISE.V.b(idxLB);
    end
end

if rompc
    [proj] = correctionProjection(ROM.H, Zbar, opt);
end

% Containers for data
z = zeros(o, T);
u = zeros(m, T);
zbar = zeros(o, T);
ubar = zeros(m, T);
tsolver = zeros(1, T);
if opt.save_xf
    xfdata = zeros(nf, T);
else
    xfdata = [];
end

% Initial conditions
xf = X0.xf;
xhat = X0.xhat;
xbar = X0.xbar;
if continuous
    xbar_k = X0.xbar;
end

% Simulate
fprintf('Percent complete:      ');
prc_cmplt = 0;
J = 0;
i = update;
for t = 1:T
    % Specify some noise terms
    wf = 0;
    v = 0;
    if isfield(opt.NOISE, 'W')
        wf = wfLB + rand(mw,1).*(wfUB - wfLB);
    elseif isfield(opt.NOISE, 'wf')
        wf = opt.NOISE.wf(:,t);
    end
    if isfield(opt.NOISE, 'V')
        v = vLB + rand(p,1).*(vUB - vLB);
    elseif isfield(opt.NOISE, 'v')
        wf = opt.NOISE.v(:,t);
    end
    
    % Get measurement and outputs
    y = FOM.Cf*xf + v;
    z(:,t) = FOM.Hf*xf;
    zbar(:,t) = ROM.H*xbar;
    if opt.save_xf
        xfdata(:,t) = xf;
    end
       
    % Solve the optimal control problem
    if rompc && precompiled
        if ~continuous || (continuous && i == update)
            if continuous
                xbar = xbar_k; % reset
            end
            [sol, ~, ~, ~, ~, solver] = ROM_CTRL(xbar);
            ubar(:,t) = sol{1};
            u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            
            % Correct for numerical errors if needed
            xbar_k = sol{2};
            if ~Zbar.contains(ROM.H*xbar_k)
                xbar_k = proj(xbar_k);
            end
            i = 1;
        elseif continuous && i ~= update
            ubar(:,t) = sol{1}; % from previous solve (ZOH)
            u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            i = i + 1;
            solver.solvertime = 0;
        end
        
        % Update state
        if continuous
            xbar = Ad*xbar + Bd*ubar(:,t);
        else
            xbar = xbar_k;
        end
        
    elseif rompc && ~precompiled
        if ~continuous || (continuous && i == update)
            if continuous
                xbar = xbar_k; % reset
            end
            [xopt, uopt, solver] = solveOptimalControl(ROM_CTRL, xbar, opt);
            ubar(:,t) = uopt(:,1);
            u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            
            % Correct for numerical errors if needed
            xbar_k = xopt(:,2);
            if ~Zbar.contains(ROM.H*xbar_k)
                xbar_k = proj(xbar_k);
            end
            i = 1;
        elseif continuous
            ubar(:,t) = uopt(:,1);
            u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            i = i + 1;
            solver.solvertime = 0;
        end
        
        % Update state
        if continuous
            xbar = Ad*xbar + Bd*ubar(:,t);
        else
            xbar = xbar_k;
        end
        
    else
        t0 = tic;
        ubar(:,t) = ROM_CTRL*xbar;
        solver.solvertime = toc(t0);
        u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
        if continuous
            xbar = Ad*xbar + Bd*ubar(:,t);
        else
            xbar = ROM.A*xbar + ROM.B*ubar(:,t);
        end
        solver.problem = 0;
    end
    
    % Check that solver succeeded
    if solver.problem ~= 0
        disp(solver);
        error('Solver failed.');
    end
    tsolver(t) = solver.solvertime;
    
    % Update the reduced order state estimate
    if continuous
        xhat = Ade*xhat + Bde*u(:,t) + Bwe*y;
    else
        xhat = ROM.A*xhat + ROM.B*u(:,t) + CTRL.L*(y - ROM.C*xhat);  
    end
    
    % Update the full order (true) system dynamics
    if continuous
        [xfsim] = simulateViaBackEuler(FOM.Af, FOM.Bf, xf, opt.sim_dt, u(:,t), zeros(nf,1), M1, M2);
        xf = xfsim(:,2);
    else
        xf = FOM.Af*xf + FOM.Bf*u(:,t) + FOM.Bfw*wf;
    end
    
    % Update cost
    if continuous
        J = J + opt.sim_dt*(xf'*FOM.Qf*xf + u(:,t)'*FOM.Rf*u(:,t));
    else
        J = J + xf'*FOM.Qf*xf + u(:,t)'*FOM.Rf*u(:,t);
    end
    
    if round(100*(t/T)) > prc_cmplt
        prc_cmplt = round(100*(t/T));
        fprintf('\b\b\b\b\b%3.0f%% ', prc_cmplt);
    end
end
fprintf('\n');

DATA.u = u;
DATA.z = z;
DATA.ubar = ubar;
DATA.zbar = zbar;
DATA.xf = xfdata;
DATA.tsolver = tsolver;
DATA.T = T;
DATA.X0 = X0;
DATA.opt = opt;
DATA.J = J;
if continuous
    DATA.t = 0:opt.sim_dt:T;
end

end

