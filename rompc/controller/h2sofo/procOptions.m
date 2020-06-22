function [pars, options] = procOptions(pars, options)

pars.penaltyConstraint = 0;

%FIRST TRANSLATE OPTIONS FROM NEW NAMES TO OLD NAMES

%  weighNormK --> penalty
%  augmentHinf --> barrier
%  rho.init --> penaltyConstraintInit
%  rho.max  --> penaltyConstraintMax
%  rho.multiply --> penaltyConstraintMultiply

if isfield(options,'weightNormK')
    options.penalty = options.weightNormK;
    options = rmfield(options,'weightNormK');
end

if isfield(options,'augmentHinf')
    options.barrier = options.augmentHinf;
    options = rmfield(options,'augmentHinf');
end

if isfield(options, 'rho')
    rho = options.rho;
    options.penaltyConstraintInit = rho.init;
    options.penaltyConstraintMultiply = rho.multiply;
    options.penaltyConstraintMax = rho.max;
    options = rmfield(options, 'rho');
end

% Number of optimization variables
pars.nvar = size(options.structure.V, 2);
if pars.nvar == 0
    % Shouldn't get here because an error would be thrown in procStructure
    error('h2sofo: no decision variables, all are fixed or there is a unique solution')
end

pars.structure = options.structure;  %already verified and set up in h2sofo.m
 
if ~isfield(options,'fast')
    options.fast = 0;
else
    if options.fast ~= 0 && options.fast ~= 1
        error('h2sofo: options.fast must be either 0 or 1')
    end
end
    
if ~isfield(options, 'penalty') 
    options.penalty = 0;
elseif isempty(options.penalty)
    options.penalty = 0;
elseif ~(isposreal(options.penalty) || options.penalty == 0)
    error('h2sofo: input "options.penalty" must be a nonnegative real scalar')
end
    
if ~isfield(options, 'barrier') 
    options.barrier = 0;
elseif isempty(options.barrier)
    options.barrier = 0;
elseif ~(isposreal(options.barrier) || options.barrier == 0)
    error('h2sofo: input "options.barrier" must be a nonnegative real scalar')
end
if ~isfield(options, 'cpumax') % quit when cpu time in seconds exceeds this
    options.cpumax = inf;
elseif ~isposreal(options.cpumax)
    error('h2sofo: input "options.cpumax" must be a positive scalar')
end
if ~isfield(options, 'normtol')
    options.normtol = 1e-3; % larger default than HANSO default
elseif ~isposreal(options.normtol)
    error('h2sofo: input "options.normtol" must be a positive scalar')
end 
if ~isfield(options, 'evaldist') 
    options.evaldist = 1e-3; 
elseif ~isposreal(options.evaldist)
    error('h2sofo: input "options.evaldist" must be a positive scalar')
end
if ~isfield(options, 'nrand') % number of starting points besides init
    options.nrand = 3;
elseif ~isposint(options.nrand)
    % nrand is 0 is not allowed: if init came from output of a hifoo run,
    % perhaps a different order, objective is likely not to be smooth
    % there and BFGS may fail immediately: some perturbation is needed
    error('h2sofo: input "options.nrand" must be a positive integer')
end
if ~isfield(options,'penaltyConstraintInit')
    options.penaltyConstraintInit = 100;
end
if ~isfield(options,'penaltyConstraintMultiply')
    options.penaltyConstraintMultiply = 100;
end

if ~isfield(options,'penaltyConstraintMax')
    options.penaltyConstraintMax = 1e10;
end

if ~isfield(options,'maxit')
    options.maxit = 1000;
elseif ~isposint(options.maxit)
    error('h2sofo: options.maxit must be a positive integer')
end     

pars.barrier = options.barrier;

end
