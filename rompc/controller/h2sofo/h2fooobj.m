function [f, g] = h2fooobj(k, pars)

%Check whether we are just stabilizing or allowing each plant to
%have its own objective function
if pars.objective == '+'
    objective = 's';
else
    objective = pars.objective;
end

if strcmp(objective, 's') % spectral abscissa
    [Acl, ~, ~, ~] = getABCDclosedloop(pars.plantinfo, k, pars.structure);
    [f, gs] = spectralabsc(Acl, pars.continuous);
    g = chainK(pars.plantinfo, gs, pars.structure);
elseif strcmp(objective, 'h2') % H2 norm
    [f, g] = h2(pars.plantinfo, k, pars.structure, pars.continuous);
else
    fprintf('h2sofo: invalid objective\n')
    f = nan;
    g = nan*ones(pars.nvar,1);
end

% the penalty term is to prevent the norm of the controller from blowing up
penalty = pars.penalty; % set by options.weightNormK I think
if penalty > 0
    knorm = norm(k);
    if knorm > 0
        f = f + penalty*knorm;
        if nargout >= 2
            g = g + penalty*k/knorm; % gradient of the 2-norm
        end
    else
        d = randn(pars.nvar);  % random subgradient in this case
        if nargout >= 2
            g = g + penalty*d/dnorm;
        end
    end
end

end