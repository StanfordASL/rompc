function [K, f, loc, options] = h2sofomain(pars, Kinit, obj, options)

pars.fgname = 'h2fooobj';
cpufinish = cputime + options.cpumax;
prtlevel = options.prtlevel;

if prtlevel > 0
    fprintf('h2sofo: number of design variables is %d\n', pars.nvar)
end
options.prtlevel = max(0, prtlevel-1); % to reduce output from HANSO

% Check out what software is available to compute H2 norm
if exist('lyap') ~= 2
    error('h2sofo: lyap function from the Control System Toolbox must be installed when input "obj" is "t"');
end

% Check whether HANSO is available: must be version 2.0 or higher which does not need quadprog   
if exist('hanso') == 2
    if prtlevel > 0
        fprintf('h2sofo: using "hanso" for optimization (should be version 2.0 or higher)\n')
    end
else
    error('h2sofo: "hanso (version 2.0)" must be installed from www.cs.nyu.edu/overton/software/hanso/hanso2_0')
end

% Initialization
if isempty(Kinit)
    % No initial point provided, starting points are randomly generated
    k0 = randn(pars.nvar, options.nrand);
    if prtlevel > 0
        fprintf('h2sofo: no initial point provided, using %d randomly generated starting points\n', options.nrand)
    end
else
    % Vectorizes Kinit
    kinit = K2k(Kinit, options.structure);
    
    % starting points include progressively bigger perturbations of xinit
    scalepert = norm(kinit)*logspace(-3, -1, options.nrand); % spaced from 0.001 to 0.1
    xinitpert = randn(pars.nvar, options.nrand)*diag(scalepert);
    k0 = [kinit, kinit*ones(1,options.nrand) + xinitpert];
    if prtlevel > 0
        fprintf('h2sofo: supplementing initial point with %d randomly generated perturbation(s)\n', options.nrand)
    end
end
oldopts = options;

if strcmp(obj, 'h2') || strcmp(obj, '+') % for H2 optimization you want to stabilize first as well
    % for these objectives, first step is to get stable
    pars.objective = '+';   % minimize spectral abscissa
    pars.penalty = 0;       % without any penalty term (barrier irrelevant)
    
    if pars.continuous
        fstable = 0;
    else
        fstable = 1;
    end
    
    options.fvalquit = fstable; % stopping when get stable
    options.maxit = 1000;   % run a long time if necessary
    if prtlevel > 0
        fprintf('h2sofo: searching for stabilizing controller\n')
        if options.nrand > 1
            fprintf('h2sofo: if stabilization time is excessive, reduce options.nrand\n')
        end
    end
    
    % use BFGS alone at first, as in most cases this will be sufficient
    % call it one starting point at a time, because when objective is '+',
    % want to stop when a stable point is found, while for objectives
    % 'h2' want to find several stable points for initializing 
    % optimization
    for i = 1:size(k0, 2)
        options.x0 = k0(:,i);
        % might waste too much of the allotted cputime stabilizing many
        % start points but complicated to avoid, user can adjust options
        options.cpumax = cpufinish - cputime;
        [kB(:,i), fB(i), gB(:,i)] = bfgs(pars, options);
        if cputime > cpufinish
            if prtlevel > 0
                fprintf('h2sofo: quit stabilizing since CPU time limit exceeded\n')
            end
            break % not return
        end
        if strcmp(obj, '+') && fB(i) < fstable % only want one stable point
            break % not return
        end
    end
    
    % if BFGS did not find a stable point, try gradient sampling,
    % initializing with the best half of the points found by BFGS,
    % starting with the lowest
    % (no point using local bundle, not trying to verify optimality)
    if cputime < cpufinish && min(fB) >= fstable
        if options.fast == 1
            if prtlevel > 1
                fprintf('h2sofo: skipping gradient sampling\n');
            end
        else
            if prtlevel > 1
                fprintf('h2sofo: BFGS did not find any stabilizing controllers, trying gradient sampling\n')
            end
            [~, indx] = sort(fB);
            indx = indx(1:ceil(length(fB)/2));
            kBsort = kB(:, indx);
            options.maxit = 100; % gradient sampling is more expensive than BFGS
            for i = 1:length(indx)
                options.k0 = kBsort(:, i);
                options.cpumax = cpufinish - cputime;
                [xGS, fGS, gGS] = gradsamp(pars, options);
                if cputime > cpufinish
                    if prtlevel > 0
                        fprintf('h2sofo: quit stabilizing since CPU time limit exceeded\n')
                    end
                    break % not return
                end
                if fGS < fstable % settle for one stable point
                    kB = xGS; fB = fGS; gB = gGS;
                    break % not return
                end
            end
        end
    end
    stabindx = find(fB < fstable);
    nstabpts = length(stabindx);
    if nstabpts > 0
        k0 = kB(:, stabindx);  % at least one stable point was found
        f0 = fB(:, stabindx);
        if prtlevel > 0
            if min(fB) > fstable - 1e-8
                qualification = '(barely) ';
            else
                qualification = '';
            end
            if nstabpts == 1
                fprintf('h2sofo: found a %sstabilizing controller', qualification)
            else
                fprintf('h2sofo: found %d %sstabilizing controllers', nstabpts, qualification)
            end
            if ~strcmp(obj, '+')
                fprintf(' for initializing optimization\n')
            else
                fprintf(' , quitting\n')
                [f, idx] = min(f0);
                ksol = k0(:, idx);
                K = k2K(pars.m, pars.p, ksol, pars.structure);
                loc = [];
                return;
            end
        end
        foundstablepoint = 1;
    else
        if prtlevel > 0
            fprintf('h2sofo: could not find a stabilizing order %d controller\n', 0)
            fprintf('h2sofo: returning controller with best spectral abscissa(radius) %g instead\n', min(fB))
            fprintf('h2sofo: try specifying a new "init" or increasing "options.cpumax"\n')
        end
        foundstablepoint = 0;
    end
    if ~foundstablepoint % nothing else to do
        K = [];
        loc.dnorm = nan;
        loc.evaldist = nan;
        return
    end
end
options = oldopts;

% i here means independent or individual.  h2fooobj will look at the set objective in pars.plantinfo
pars.objective = obj;
switch obj
    case 's'
        if pars.continuous
            objectivename = 'spectral abscissa';
        else
            objectivename = 'spectral radius';
        end
    case '+'
        objectivename = 'stabilize only';
    case 'h2'
	    objectivename = 'H2 performance';	
    otherwise
        objectivename = [];
end
pars.objname = objectivename; % may be handy for display purposes
pars.penalty = options.penalty;   % penalty term on ||x||_2
if pars.penalty > 0
    penaltystring = ' (plus penalty term)';
else
    penaltystring = '';
end
pars.barrier =  options.barrier;   % barrier term: multiple of inverse of stability radius
if pars.barrier > 0
    barrierstring = ' (plus barrier term)';
else
    barrierstring = '';
end
options.fvalquit = -inf; % do not quit early
options.x0 = k0;   % multiple starting points, stable if obj is 'h', 't' or 'r'
if prtlevel > 0
    fprintf('h2sofo: optimizing %s%s%s for order %d controller\n', ...
        objectivename, penaltystring, barrierstring, 0)
    fprintf('h2sofo: if optimization time is excessive, reduce options.cpumax and/or options.nrand\n')
end
% run HANSO from all the valid starting points, returning only best result
% repeat running hanso with larger and larger penalty mutlipliers
% until we find a solution which satisfies all constraints

% pars.penaltyConstraintInit = options.penaltyConstraintInit;
% pars.penaltyConstraint = options.penaltyConstraintInit;
options.cpumax = cpufinish - cputime; % time left after stabilization
cpufinish = cputime + options.cpumax;

options.samprad = []; % turn off gradient sampling (HANSO 2.0)
% in case of HANSO 1.0 is called, turn off local bundle and gradient sampling
options.phasenum = [size(options.x0, 2) 0 0]; 
[ksol, ~, loc] = hanso(pars, options);

if options.fast == 0 
    if prtlevel > 1
        fprintf('h2sofo: attempting to improve solution via gradient sampling\n');
    end
    % RUN GRADIENT SAMPLING (via HANSO) TO IMPROVE SOLUTION FROM BFGS ABOVE
    % Only run from the best point above.
    options.k0 = ksol;
    options.cpumax = cpufinish - cputime;
    % options.phasenum = [0,1,2]; was used by HANSO 1.0 to tell it to skip BFGS.
    % not used by HANSO 2.0, which will try BFGS and then immediately
    % switch to gradient sampling when BFGS fails as it has already been run
    % But gradient sampling is not done if the number of variables > 100 
    % as it is too expensive: then nothing will happen so xsol will not change
    options = rmfield(options, 'samprad'); % turn gradient sampling back on
    options = rmfield(options,'phasenum');  % turn gradient sampling back on (hanso 1.0)
    [ksol, ~, loc] = hanso(pars,options);
end

f = h2fooobj(ksol, pars);

% Verify and issue a warning if any of the closed loop systems are
% approaching the stability boundary
checkStabilityBoundary(ksol, pars, prtlevel);

if isnan(loc.dnorm)
    loc.dnorm = zeros(size(loc.dnorm));
end

if prtlevel > 0       
    fprintf('h2sofo: best order %d controller found has %s%s%s %g\n', ...
        0, objectivename, penaltystring, barrierstring, f)

    if exist('loc','var')
        fprintf(' with local optimality measure: dnorm = %5.1e, evaldist = %5.1e\n',...
            loc.dnorm, loc.evaldist')
    end
end

% Return the best point found, translated into controller format
K = k2K(pars.m, pars.p, ksol, pars.structure);

end
