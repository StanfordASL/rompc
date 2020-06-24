function [Xbar] = computeXbar(A, B, H, Z, U, k, idx, opt)
%[Xbar] = computeXbar(A, B, H, Z, U, k, idx, opt)
%
%Computes a bounding hyperrectangular polytope for the ROM state xbar at
%time t-k based on the assumption that H*xbar \in Z and ubar \in U in 
%interval [t-k, t] (i.e. past k time steps).
%
%The system A,B,H is a discrete time system, and therefore k and t are time
%steps and not seconds
%
% Inputs:
%   A: discrete time reduced order system matrix
%   B: discrete time reduced order control matrix
%   H: reduced order performance matrix
%   Z: performance variable constraints (Polyhedron)
%   U: control constraints (Polyhedron)
%   k: time horizon to consider
%   idx: value in range [1, k] for which state to bound in [t-k,k],
%        Lorenzetti (2020) paper uses t = 1 as default for assumptions.
%   opt: optional arguments, including:
%       - solver: specifies the solver to use
%       - xbar_datapath: path to file where precomputed matrices are
%                             stored to or loaded from
%
% Outputs:
%   Xbar: polyhedron outer bound for xbar

if ~isfield(opt, 'solver')
    opt.solver = 'cplex';
end

% Check whether the matrices that are computed should be saved or loaded
if isfield(opt, 'xbar_datapath')

    % Load data from file
    if exist(opt.xbar_datapath, 'file')
        fprintf('******************************\n');
        fprintf('Loading Xbar from %s.\nDelete this file for recomputation.\n', opt.xbar_datapath);
        fprintf('******************************\n\n');
        load(opt.xbar_datapath, 'Xbar', 'k', 'idx');
        fprintf('Using data values of k = %d and t = %d.\n', k, idx);
        return;

    % Compute data matrices and then save to file
    elseif ~exist(opt.xbar_datapath, 'file')
        fprintf('Will save data to %s\n', opt.xbar_datapath);
        savedata = true;
    end

% Compute data matrices and do not save
else
    savedata = false;
end

fprintf('Computing Xbar with k = %d and idx = %d.\n', k, idx);
n = size(A,1);
m = size(B,2);
if k < n
    fprintf('k must be >= %d to ensure boundedness.\n', n);
    Xbar = 0;
    return
end

% Checking observability of H on xbar
observable = rank(obsv(A, H)) == n;
if observable
    fprintf('The pair (A,H) is observable.\n');
else
    fprintf('The pair (A,H) may not be observable (by rank test).\n');
end

% Define decision variables
xbar = {};
for i = 1:k+1 % 0 to k
    xbar{i} = sdpvar(n,1);
end

ubar = {};
for i = 1:k % 0 to k-1
    ubar{i} = sdpvar(m,1); 
end

% Constraints
constraints = [];
    
% Dynamics constraints
for i = 1:k % 0 to k-1
    constraints = [constraints, xbar{i+1} - A*xbar{i} - B*ubar{i} == 0]; % dynamics of ROM
end

% Performance and control constraints
for i = 1:k+1 % 0 to k
    constraints = [constraints, Z.A*H*xbar{i} <= Z.b]; % reduced performance constraints
end

for i = 1:k % 0 to k-1
    constraints = [constraints, U.A*ubar{i} <= U.b]; % control constraints
end

% Just do box constraint
P = [eye(n); -eye(n)];

% Solve LPs to get p^Tx <= b
ops = sdpsettings('verbose', 0, 'solver', opt.solver, 'savesolveroutput', 1);
n2 = size(P,1);
b = zeros(n2, 1);
for i = 1:n2
    fprintf(1, 'Computing %d of %d \n', i, n2);
    J = -P(i,:)*xbar{idx};
    diagnostics = optimize(constraints, J, ops);
    b(i) = -value(J);
   
    % Cplex results
    if diagnostics.problem ~= 0
        disp('Solver failed.\n');
        disp(diagnostics.info)
        b(i) = inf;
    end
end

Xbar = Polyhedron(P, b);

% Save the data
if savedata
    fprintf('Saving data to %s.\n\n', opt.xbar_datapath);
    save(opt.xbar_datapath, 'Xbar','k','idx');
end
end

