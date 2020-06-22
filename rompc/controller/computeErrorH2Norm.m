function [f] = computeErrorH2Norm(FOM, ROM, CTRL, opts)
%[f] = computeErrorH2Norm(FOM, ROM, CTRL, opts)
%
%Compute the H2 norm of the error system.
%

fprintf('Computing H2 norm of the error system.\n');
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

nf = size(FOM.Af,1);
n = size(ROM.A,1);
m = size(FOM.Bf,2);
p = size(FOM.Cf,1);
o = size(FOM.Hf,1);
Pp = eye(nf) - ROM.V*inv(ROM.W'*ROM.V)*ROM.W';
Ae = [FOM.Af, FOM.Bf*CTRL.K; 
      CTRL.L*FOM.Cf, ROM.A + ROM.B*CTRL.K - CTRL.L*ROM.C];
Be = [Pp*FOM.Af*ROM.V, Pp*FOM.Bf; 
      zeros(n, n), zeros(n, m)];
  
if opts.with_noise && isfield(FOM, 'Bfw') && ~isempty(FOM.Bfw)
    fprintf('Including external noise terms in H2 norm calculation.\n');
    mw = size(FOM.Bfw, 2);
    Be = [Be, [FOM.Bfw, zeros(nf, p); zeros(n, mw), CTRL.L]];
end

He = [FOM.Hf, zeros(o, n);
      zeros(m, nf), CTRL.K];
Dzu = zeros(size(He,1), size(Be,2));
if opts.continuous
    sys = ss(full(Ae), full(Be), He, Dzu);
else
    sys = ss(full(Ae), full(Be), He, Dzu, -1);
end
f = norm(sys, 2);

end

