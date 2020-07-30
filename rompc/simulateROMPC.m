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
%       - rom_only: true to only simulate the ROM, default = false
%
% Returns:
%   DATA: struct containing simulation data as well as options and
%         controller query times

if ~isfield(opt, 'rom_only')
    opt.rom_only = false;
end

if opt.rom_only || ~isfield(opt, 'save_xf')
    opt.save_xf = false;
end

if ~opt.rom_only
    nf = size(FOM.Af, 1);
end
m = size(ROM.B, 2);
p = size(ROM.C, 1);
o = size(ROM.H, 1);

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
    if ~opt.rom_only
        S = speye(nf) - opt.sim_dt*FOM.Af;
        setup.type = 'nofill';
        [M1, M2] = ilu(sparse(S), setup);
    end
    
elseif isfield(opt, 'discrete')
    continuous = ~opt.discrete;
    update = -1;
else
    error('Need to specify if the problem is discrete or continuous time.');
end

% Check what kind of controller was used for ROM
if isstruct(ROM_CTRL)
    % Open loop control
    if isfield(ROM_CTRL, 'ubar') && isfield(ROM_CTRL, 'xbar')
        ctrl_type = 'openloop';
        X0.xbar = ROM_CTRL.xbar(:,1);
        if size(ROM_CTRL.ubar, 2) == T-1
            ROM_CTRL.ubar = [ROM_CTRL.ubar, ROM_CTRL.ubar(:,end)];
        end
        if ~isfield(ROM_CTRL, 't') && continuous
            error('Open loop needs time field for continuous time problems.');
        end

    % ROMPC without precompiled optimizer
    else
        ctrl_type = 'rompc';
    end
elseif isa(ROM_CTRL, 'optimizer')
    ctrl_type = 'rompc_precompiled';
elseif ismatrix(ROM_CTRL)
    ctrl_type = 'linearfeedback';
end

if ~opt.rom_only
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
else
    opt.NOISE = struct();
end

if strcmp(ctrl_type, 'rompc') || strcmp(ctrl_type, 'rompc_precompiled')
    [proj] = correctionProjection(ROM.H, Zbar, opt);
end

% Containers for data
if ~opt.rom_only
    z = zeros(o, T);
    u = zeros(m, T);
end
zbar = zeros(o, T);
ubar = zeros(m, T);
tsolver = zeros(1, T);
if opt.save_xf
    xfdata = zeros(nf, T);
else
    xfdata = [];
end

% Initial conditions
if ~opt.rom_only
    xf = X0.xf;
    xhat = X0.xhat;
end
xbar = X0.xbar;
if continuous
    xbar_k = X0.xbar;
end

% Simulate
fprintf('Percent complete:      ');
prc_cmplt = 0;
J = 0;
Jr = 0;
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
    if ~opt.rom_only
        y = FOM.Cf*xf + v;
        z(:,t) = FOM.Hf*xf;
        if opt.save_xf
            xfdata(:,t) = xf;
        end
    end
    zbar(:,t) = ROM.H*xbar;
     
    % Solve the optimal control problem
    if strcmp(ctrl_type, 'rompc_precompiled')
        if ~continuous || (continuous && i == update)
            if continuous
                xbar = xbar_k; % reset
            end
            [sol, ~, ~, ~, ~, solver] = ROM_CTRL(xbar);
            ubar(:,t) = sol{1};
            if ~opt.rom_only
                u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            end
            
            % Correct for numerical errors if needed
            xbar_k = sol{2};
            if ~Zbar.contains(ROM.H*xbar_k)
                xbar_k = proj(xbar_k);
                disp(xbar_k);
            end
            i = 1;
        elseif continuous && i ~= update
            ubar(:,t) = sol{1}; % from previous solve (ZOH)
            if ~opt.rom_only
                u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            end
            i = i + 1;
            solver.solvertime = 0;
        end
        
        % Update state
        if continuous
            xbar = Ad*xbar + Bd*ubar(:,t);
        else
            xbar = xbar_k;
        end
        
    elseif strcmp(ctrl_type, 'rompc')
        if ~continuous || (continuous && i == update)
            if continuous
                xbar = xbar_k; % reset
            end
            [xopt, uopt, solver] = solveOptimalControl(ROM_CTRL, xbar, opt);
            ubar(:,t) = uopt(:,1);
            if ~opt.rom_only
                u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            end
            
            % Correct for numerical errors if needed
            xbar_k = xopt(:,2);
            if ~Zbar.contains(ROM.H*xbar_k)
                xbar_k = proj(xbar_k);
            end
            i = 1;
        elseif continuous
            ubar(:,t) = uopt(:,1);
            if ~opt.rom_only
                u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            end
            i = i + 1;
            solver.solvertime = 0;
        end
        
        % Update state
        if continuous
            xbar = Ad*xbar + Bd*ubar(:,t);
        else
            xbar = xbar_k;
        end
        
    elseif strcmp(ctrl_type, 'linearfeedback')
        t0 = tic;
        ubar(:,t) = ROM_CTRL*xbar;
        solver.solvertime = toc(t0);
        if ~opt.rom_only
            u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
        end
        if continuous
            xbar = Ad*xbar + Bd*ubar(:,t);
        else
            xbar = ROM.A*xbar + ROM.B*ubar(:,t);
        end
        solver.problem = 0;

    elseif strcmp(ctrl_type, 'openloop')
        if continuous
            time = (t-1)*opt.sim_dt;
            timenext = t*opt.sim_dt;
            if timenext > ROM_CTRL.t(end)
                break;
            end
            ubar(:,t) = interp1(ROM_CTRL.t, ROM_CTRL.ubar', time);
            if ~opt.rom_only
                u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            end
            xbar = interp1(ROM_CTRL.t, ROM_CTRL.xbar', timenext);
            xbar = xbar';
        else
            if t >= size(ROM_CTRL.ubar, 2)
                break;
            end
            ubar(:,t) = ROM_CTRL.ubar(:,t);
            if ~opt.rom_only
                u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
            end
            xbar = ROM_CTRL.xbar(:,t+1);
            xbar = xbar';
        end       
        solver.solvertime = 0;
        solver.problem = 0;
    end
    
    % Check that solver succeeded
    if solver.problem ~= 0
        disp(solver);
        error('Solver failed.');
    end
    tsolver(t) = solver.solvertime;
    
    if ~opt.rom_only
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
    end
    
    % Update reduced order cost
    if continuous
        Jr = Jr + opt.sim_dt*(xbar'*ROM.Q*xbar + ubar(:,t)'*ROM.R*ubar(:,t));
    else
        Jr = Jr + xbar'*ROM.Q*xbar + ubar(:,t)'*ROM.R*ubar(:,t);
    end
    
    if round(100*(t/T)) > prc_cmplt
        prc_cmplt = round(100*(t/T));
        fprintf('\b\b\b\b\b%3.0f%% ', prc_cmplt);
    end
end
fprintf('\n');

if ~opt.rom_only
    DATA.u = u;
    DATA.z = z;
    DATA.xf = xfdata;
end
DATA.ubar = ubar;
DATA.zbar = zbar;
DATA.tsolver = tsolver;
DATA.T = T;
DATA.X0 = X0;
DATA.opt = opt;
DATA.J = J;
DATA.Jr = Jr;
if continuous
    DATA.t = opt.sim_dt*(0:T-1);
end

end

