function [DATA] = simulateROMPC(FOM, ROM, CTRL, ROMPC, Zbar, T, X0, opt)
%[DATA] = simulateROMPC(FOM, ROM, CTRL, ROMPC, Zbar, T, X0, opt)
%
%Simulate the FOM with ROMPC control.
%
% Inputs:
%   FOM: struct containing FOM matrices (Af, Bf, Bfw, Cf, Hf)
%   ROM: struct containing ROM matrices (A, B, Bw, C, H, V, W)
%   CTRL: struct containing gains K and L
%   ROMPC: controller for ROM, either is the precompiled ROMPC problem, a
%          structure containing the ROMPC problem information, or a gain
%          matrix for state feedback control
%   Zbar: tightened constraints used in ROMPC
%   T: time horizon
%   X0: struct containing xf, xhat, xbar initial conditions
%   opt: optional arguments, including
%       - NOISE: struct containing W or V Polyhedrons for random bounded
%                disturbances, or wf/v matrices of size mw/p by T to
%                directly specify the noise terms at each time step
%
% Returns:
%   DATA: struct containing simulation data as well as options and
%         controller query times

nf = size(FOM.Af, 1);
m = size(FOM.Bf, 2);
p = size(FOM.Cf, 1);
o = size(FOM.Hf, 1);

% Check what kind of controller was used for ROM
if isstruct(ROMPC)
    precompiled = false;
    rompc = true;
elseif strcmp(class(ROMPC), 'optimizer')
    precompiled = true;
    rompc = true;
elseif ismatrix(ROMPC)
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

% Containers for data
z = zeros(o, T);
u = zeros(m, T);
zbar = zeros(o, T);
ubar = zeros(m, T);
tsolver = zeros(1,T);

% Initial conditions
xf = X0.xf;
xhat = X0.xhat;
xbar = X0.xbar;

% Simulate
fprintf('Percent complete:      ');
prc_cmplt = 0;
J = 0;
[proj] = correctionProjection(ROM.H, Zbar, opt);
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
       
    % Solve the optimal control problem
    if rompc && precompiled 
        [sol, ~, ~, ~, ~, solver] = ROMPC(xbar);
        ubar(:,t) = sol{1};
        u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
        xbar = sol{2};
        
        % Correct for numerical errors if needed
        if ~Zbar.contains(ROM.H*xbar)
            xbar = proj(xbar);
        end
    elseif rompc && ~precompiled
        [xopt, uopt, solver] = solveOptimalControl(ROMPC, xbar, opt);
        ubar(:,t) = uopt(:,1);
        u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
        xbar = xopt(:,2);
        
        % Correct for numerical errors if needed
        if ~Zbar.contains(ROM.H*xbar)
            xbar = proj(xbar);
        end
    else
        t0 = tic;
        ubar(:,t) = ROMPC*xbar;
        solver.solvertime = toc(t0);
        u(:,t) = ubar(:,t) + CTRL.K*(xhat - xbar);
        xbar = ROM.A*xbar + ROM.B*ubar(:,t);
        solver.problem = 0;
    end
    
    % Check that solver succeeded
    if solver.problem ~= 0
        disp(solver);
        error('Solver failed.');
    end
    tsolver(t) = solver.solvertime;
    
    % Update the reduced order state estimate
    xhat = ROM.A*xhat + ROM.B*u(:,t) + CTRL.L*(y - ROM.C*xhat);  
    
    % Update the full order (true) system dynamics
    xf = FOM.Af*xf + FOM.Bf*u(:,t) + FOM.Bfw*wf;
    J = J + xf'*FOM.Qf*xf + u(:,t)'*FOM.Rf*u(:,t);
    
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
DATA.tsolver = tsolver;
DATA.T = T;
DATA.X0 = X0;
DATA.opt = opt;
DATA.J = J;
end

