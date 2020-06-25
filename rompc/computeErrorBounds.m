function [EBOUND] = computeErrorBounds(FOM, ROM, CTRL, Z, U, Wnoise, Vnoise, tau, opt)
%[EBOUND] = computeErrorBounds(FOM, ROM, CTRL, Z, U, Wnoise, Vnoise, tau, opt)
%
%Computes error bound parameters for ROMPC.
%
% Inputs:
%   FOM: struct containing FOM matrices (Af, Bf, Bfw, Cf, Hf)
%   ROM: struct containing ROM matrices (A, B, Bw, C, H, V, W)
%   CTRL: struct containing gains K and L
%   Z: performance variable constraints (Polyhedron)
%   U: control constraints (Polyhedron)
%   Wnoise: process noise bound (Polyhedron), or 0 if no noise
%   Vnoise: measurement noise bound (Polyhedron), or 0 if no noise
%   tau: backward time horizon
%   opt: options, including:
%       - solver: string, solver to use e.g. 'cplex' or 'mosek'
%       - G_method:   'Lyap', 'LMI', 'GP', 'BatchGP', or 'None'
%       - G_0:        PD matrix to warm start G_method for 'BatchGP' (diagonal), or 
%                     to use if 'None', size nf + n x nf + n
%       - Xbar_tau:   alternate time horizon for computing Xbar
%       - Xbar_t:     default t = 1, can be up to t = Xbar_tau
%       - Cr_method:   'upperbound', 'hyperrect', 'default'
%       - Cw_method:   'upperbound', 'hyperrect', 'default'
%       - continuous: true to specify the FOM, ROM, CTRL are for a continuous
%                     time system
%       - normbound_datapath: path for where to save or load normbound
%                             data from
%       - errormats_datapath: path for where to save or load error matrix
%                             data used in computing D2
%       - xbar_datapath:      path for where to save or load xbar
%       - D2_only:    true to only compute D2 error terms
%       - no_FOM:     true if don't have the FOM but only data matrices
%                     saved in errormats_datapath
%       - dt:             For continuous time systems, specify the time discretization
%                         to use in quadrature scheme
%       - int_method:     For continuous time systems, specify the quadrature scheme
%                         options include 'trapezoid' or 'leftRiemann'
%       - mat_exp_method: For continuous time systems, specify method for approximating
%                         the matrix exponential terms, either 'expm' for direct computation
%			              'expmv' for approximation of mat exp action, or 'implicitEuler'
%			              to use an implicit Euler integration scheme
%
% Outputs:
%   EBOUND: struct containing relevant parameters for error bound

% Extract options
if isfield(opt, 'continuous')
    continuous = opt.continuous;
elseif isfield(opt, 'discrete')
    continuous = ~opt.discrete;
else
    fprintf('Need to specify if the problem is discrete or continuous time.\n');
    return;
end

if isfield(opt, 'D2_only')
    D2_only = opt.D2_only;
else
    D2_only = false;
end

if isfield(opt, 'no_FOM')
    no_FOM = opt.no_FOM;
    if no_FOM
        D2_only = true;
    end
else
    no_FOM = false;
end

if continuous
    if isfield(opt, 'dt')
        dt = opt.dt;
        fprintf('Considering a continuous time problem.\n');
    else
        fprintf('Continuous time problem, but no dt specified, using dt = 0.01\n');
        dt = 0.01;
    end
else
    fprintf('Considering a discrete time problem.\n');
end

if isfield(opt, 'G_method')
    G_method = opt.G_method;
else
    G_method = 'None';
    fprintf('Using G_method: None\n');
end

if isfield(opt, 'G_0')
    G_0 = opt.G_0;
else
    G_0 = [];
end

if isfield(opt, 'Cr_method')
    Cr_method = opt.Cr_method;
else
    Cr_method = 'none';
end

if isfield(opt, 'Cw_method')
    Cw_method = opt.Cw_method;
else
    Cw_method = 'none';
end

if isfield(opt, 'Xbar_tau')
    Xbar_tau = opt.Xbar_tau;
else
    Xbar_tau = tau;
    fprintf('Xbar_tau not specified, using: %d\n', Xbar_tau);
end

if isfield(opt, 'Xbar_t')
    Xbar_t = opt.Xbar_t;
else
    if continuous
        Xbar_t = 0;
    else
        Xbar_t = 1;
    end
end

if continuous
    T = max(floor(Xbar_tau/dt), size(ROM.A,1));
    Xbar_idx = floor(Xbar_t/dt) + 1;
    if Xbar_idx == T + 1
        Xbar_idx = T;
    end
else
    T = max(Xbar_tau, size(ROM.A,1));
    Xbar_idx = Xbar_t;
end

if Xbar_idx < 1 || Xbar_idx > T
    fprintf('Xbar_t must be in range [1,%d].\n', T);
    return;
end

% Compute Xbar 
if continuous
    [A, B, ~] = zoh(dt, ROM.A, ROM.B, 0);
else
    A = ROM.A; B = ROM.B;
end
[Xbar] = computeXbar(A, B, ROM.H, Z, U, T, Xbar_idx, opt);


% Error system
if no_FOM
    ERROR = [];
else
    [ERROR] = buildErrorSystem(FOM, ROM, CTRL, Z, U, Wnoise, Vnoise);
end

% Compute D2
if continuous
    [Dz_2, Du_2] =  computeInputEffectBoundContinuous(ROM, ERROR, Z, U, Xbar, Wnoise, Vnoise, tau, opt);
else
    [Dz_2, Du_2] =  computeInputEffectBoundDiscrete(ROM, ERROR, Z, U, Xbar, Wnoise, Vnoise, tau, opt);
end

% Compute D2 only, this also defaulted to when no_FOM = true
if D2_only
    if continuous
        EBOUND = struct('Dz_2', Dz_2, 'Du_2', Du_2, 'Xbar', Xbar, 'tau', tau, 'dt', dt);
    else
        EBOUND = struct('Dz_2', Dz_2, 'Du_2', Du_2, 'Xbar', Xbar, 'tau', tau);
    end

    return;
end

% Compute norm weighting matrix G and M, gamma
if continuous
    [beta, alpha, G, sqrtG] = computeScalarDecayBoundContinuous(ERROR, G_method, G_0, opt);
    fprintf('Using beta = %.3f, alpha = %.3f.\n', beta, alpha);
else
    [M, gamma, G, sqrtG] = computeScalarDecayBoundDiscrete(ERROR, G_method, G_0, opt);
    fprintf('Using M = %.3f, gamma = %.3f.\n', M, gamma);
end

% Compute Cr and Cw
if strcmp(Cr_method, 'upperbound')
    fprintf('Computing upper bound on Cr using convex relaxation.\n');
    [Cr] = computeC_UpperBound(ERROR.Be, G, Xbar, U, opt);
elseif strcmp(Cr_method, 'hyperrect')
    fprintf('Computing Cr assuming Xbar and U are hyperrectangles.\n');
    [Cr] = computeC_HyperRect(ERROR.Be, G, Xbar, U);
else
    fprintf('Computing Cr by checking all vertices (not assuming hyperrectangle).\n');
    [Cr] = computeC_VertexEnum(ERROR.Be, G, Xbar, U);
end
    
if isa(Wnoise, 'Polyhedron') && isa(Vnoise, 'Polyhedron')
    if strcmp(Cw_method, 'upperbound')
        fprintf('Computing upper bound on Cw using convex relaxation.\n');
        [Cw] = computeC_UpperBound(ERROR.Ge, G, Wnoise, Vnoise, opt);
    elseif strcmp(Cw_method, 'hyperrect')
        fprintf('Computing Cw assuming Xbar and U are hyperrectangles.\n');
        [Cw] = computeC_HyperRect(ERROR.Ge, G, Wnoise, Vnoise);
    else
        fprintf('Computing Cw by checking all vertices (not assuming hyperrectangle).\n');
        [Cw] = computeC_VertexEnum(ERROR.Ge, G, Wnoise, Vnoise);
    end
else
    Cw = 0;
end

% Compute parts of D1 norm bound
if continuous
    D_1a = beta*exp(2*alpha*tau);
    D_1b = beta*exp(alpha*tau)*(Cr+Cw)/(-alpha);
else
    D_1a = M*gamma^(2*tau);
    D_1b = M*(gamma^tau)*(Cr+Cw)/(1-gamma);
end

% Compute combined other parameters
nz = size(Z.b, 1);
eta_z = zeros(nz, 1);
invsqrtG = inv(sqrtG);
for i = 1:nz
   eta_z(i) =  norm(ERROR.Ez(i,:)*invsqrtG);
end

nu = size(U.b, 1);
eta_u = zeros(nu, 1);
for i = 1:nu
    eta_u(i) = norm(ERROR.Eu(i,:)*invsqrtG);
end

% Output all data
if continuous
    EBOUND = struct('Dz_2', Dz_2, 'Du_2', Du_2, 'eta_z', eta_z, 'eta_u', eta_u, 'D_1a', D_1a, 'D_1b', D_1b, ...
            'Cr', Cr, 'Cw', Cw, 'Xbar', Xbar, 'beta', beta, 'alpha', alpha, 'G', sparse(G), 'tau', tau, 'dt', dt);
else
    EBOUND = struct('Dz_2', Dz_2, 'Du_2', Du_2, 'eta_z', eta_z, 'eta_u', eta_u, 'D_1a', D_1a, 'D_1b', D_1b, ...
            'Cr', Cr, 'Cw', Cw, 'Xbar', Xbar, 'M', M, 'gamma', gamma, 'G', sparse(G), 'tau', tau);
end

end
