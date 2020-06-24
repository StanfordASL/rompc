function [K, L] = computeControllerGains(FOM, ROM, opts)
%[K, L] = computeControllerGains(FOM, ROM, opts)
%
%Computes an approximately H2 optimal controller for the ROM tracking
%problem
%
% Inputs:
%   FOM: struct containing FOM matrices (Af, Bf, Bfw, Cf, Hf)
%   ROM: struct containing ROM matrices (A, B, Bw, C, H, V, W)
%   opts: additional options
%       - continuous (discrete): can be true or false, this is required
%       - method: can be 'h2sofo' or 'ric'
%       - with_noise: true to include FOM.Bfw terms
%       - warmstart: true or false, only used if method = 'h2sofo'
%       - Wz, weights for z, default = I
%       - Wu, weights for u, default = I
%       - regularize: add some regularization with method = 'ric', makes
%         cost H'*Wz'*Wz*H + regularize*I, default = 0
%
% Returns:
%   K: controller gains
%   L: estimator gains

n = size(ROM.A,1);
m = size(ROM.B,2);
p = size(ROM.C,1);
o = size(ROM.H,1);

if ~isfield(opts, 'continuous') && ~isfield(opts, 'discrete')
    error('Need to specify continuous or discrete system.');
elseif isfield(opts, 'continuous') && isfield(opts, 'discrete')
    if opts.continuous == opts.discrete
        error('Only specify opts.continuous or opts.discrete, not both.');
    end
elseif isfield(opts, 'discrete')
    opts.continuous = ~opts.discrete;
end

if ~isfield(opts, 'with_noise')
    opts.with_noise = false;
end

if ~isfield(opts, 'method')
    opts.method = 'ric';
end

if ~isfield(opts, 'warmstart')
    opts.warmstart = false;
end

if ~isfield(opts, 'Wz')
    opts.Wz = eye(o);
end

if ~isfield(opts, 'Wu')
    opts.Wu = eye(m);
end

if ~isfield(opts, 'regularize')
    opts.regularize = 0;
end

% Riccati based method
if strcmp(opts.method, 'ric')
    fprintf('Computing the LQG solution via model reduction.\n');
    QK = ROM.H'*opts.Wz'*opts.Wz*ROM.H + opts.regularize*eye(n);
    RK = opts.Wu'*opts.Wu;
    if isfield(ROM, 'Bw') && size(ROM.Bw,1) == n && opts.with_noise
        QL = ROM.Bw*ROM.Bw' + eye(n);
    else
        QL = eye(n);
    end
    RL = eye(p);

    if opts.continuous
        [~,~,K] = care(ROM.A, -ROM.B, QK, RK);
        [~,~,L] = care(ROM.A', ROM.C', QL, RL);
        L = L';
    else
        [~,~,K] = dare(ROM.A, -ROM.B, QK, RK);
        [~,~,L] = dare(ROM.A', ROM.C', QL, RL);
        L = L';
    end
    return;
end

% If not Riccati method, then need
nf = size(FOM.Af,1);
mw = size(FOM.Bfw,2);
Pp = eye(nf) - ROM.V*inv(ROM.W'*ROM.V)*ROM.W';
Be = [Pp*FOM.Af*ROM.V, Pp*FOM.Bf; 
      zeros(n, n), zeros(n, m)];
A_tilde = blkdiag(FOM.Af, ROM.A);
B_tilde = [zeros(nf,n), FOM.Bf;
          eye(n), ROM.B];
C_tilde = [FOM.Cf, -ROM.C;
          zeros(n, nf), eye(n)];
H_tilde = [opts.Wz*FOM.Hf, zeros(size(opts.Wz,1), n);
           zeros(size(opts.Wu,1), nf), zeros(size(opts.Wu,1), n)];
Dzu_tilde = [zeros(o, n), zeros(o, m);
             zeros(size(opts.Wu,1), n), opts.Wu];
ctrl_structure = [nan(n,p), zeros(n,n);
                  zeros(m,p), nan(m,n)];

if strcmp(opts.method, 'h2sofo')
    fprintf('Computing controller with H2SOFO.\n');
    if opts.with_noise && isfield(FOM, 'Bfw') && ~isempty(FOM.Bfw)
        fprintf('Found noise terms, including them in H2 optimization.\n');
        Bw_tilde = [Be, [FOM.Bfw, zeros(nf,p);
                         zeros(n, mw), zeros(n,p)]];
        Dyw_tilde = [zeros(p, n), zeros(p, m), zeros(p, mw), eye(p);
                     zeros(n, n), zeros(n, m), zeros(n, mw), zeros(n, p)];
    else
        Bw_tilde = Be;
        Dyw_tilde = [zeros(p, n), zeros(p, m);
                     zeros(n, n), zeros(n, m)];
    end
    
    P = struct('A', A_tilde, 'B1', Bw_tilde, 'B2', B_tilde, ...
                'C1', H_tilde, 'C2', C_tilde, ...
                'D11', zeros(size(H_tilde,1), size(Bw_tilde,2)), ...
                'D12', Dzu_tilde ,'D21', Dyw_tilde);

    bfgsoptions.fast = 1;
    bfgsoptions.prtlevel = 3;
    bfgsoptions.continuous = opts.continuous;
    bfgsoptions.structure = ctrl_structure;

    if isfield(opts, 'init')
        fprintf('Initializing controller with user provided values.\n');
        bfgsoptions.init = opts.init;
    elseif opts.warmstart
        fprintf('Using H2SOFO but initializing controller by solving an LQG problem\n');
        fprintf('via the reduced order model.\n');
        warmstart_opts = opts;
        warmstart_opts.method = 'ric';
        [Kinit, Linit] = h2optController(FOM, ROM, warmstart_opts);
        bfgsoptions.init = [Linit, zeros(n,n);
                            zeros(m,p), Kinit];
    end
    
    [K_tilde, f] = h2sofo(P, 'h2', bfgsoptions);
    fprintf('H2 norm is %.2f.\n', f);
    L = K_tilde(1:n, 1:p);
    K = K_tilde(n+1:end, p+1:end);
end

end

