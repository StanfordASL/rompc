function [beta, alpha, G, sqrtG] = computeScalarDecayBoundContinuous(ERROR, G_method, G_0, opt)
%[beta, alpha, G, sqrtG] = computeScalarDecayBoundContinuous(ERROR, G_method, G_0, opt)
%
%Computes a bound ||e^At||_G \leq beta*e^alpha*t for continuous time problems
%
% Inputs:
%   ERROR: error system struct with Ae, Be, Ge, Ez, and Eu matrices
%   G_method: method for computing a good norm weighting, either 'Lyap',
%             'LMI', 'GP', 'BatchGP', or 'None'
%   G_0: an initial guess if method is 'BatchGP' or to just use if method
%        is 'None'
%   opt: optional arguments including
%       - solver: specifies solver to use for 'LMI', 'GP', and 'BatchGP'
%       - normbound_datapath: path to save or load results
%
% Returns:
%   beta: constant in norm bound
%   alpha: exponential term constant in norm bound
%   G: PD matrix defining weighted norm x^TGx = ||x||_G^2
%   sqrtG: matrix square root of G

ne = size(ERROR.Ae, 1);

% First check if data can be loaded, or if it should be saved
if isfield(opt, 'normbound_datapath')
    if exist(opt.normbound_datapath, 'file')
        fprintf('******************************\n');
        fprintf('Loading norm transient bound alpha, beta from %s.\nDelete this file to recompute.\n', opt.normbound_datapath);
        fprintf('******************************\n\n');
        load(opt.normbound_datapath, 'beta', 'alpha', 'G', 'sqrtG');
        G = full(G);
        sqrtG = full(sqrtG);
        return;
    elseif ~exist(opt.normbound_datapath, 'file')
        fprintf('Will save norm bound data to %s\n', opt.normbound_datapath);
        savedata = true;
    end
else
    savedata = false;
end

% Compute the required quantites
if strcmp(G_method, 'log-norm')
    fprintf('Computing G, beta, alpha via logarithmic norm.\n');
    if all(size(G_0) == [ne, ne])
        fprintf('Not computing G, using G_0 as specified. Checking if G_0 is PD.\n');
        % Check to make sure G_0 is positive definite
        [~, p] = chol(G_0);
        if p ~= 0
            fprintf('G_0 is not PD.\n');
            return;
        end
        G = G_0;
    else
        fprintf('G_0 incorrect dimension, using G_0 = I instead.\n');
        G = eye(ne);
    end
    sqrtG = sqrtm(G);
    invsqrtG = inv(sqrtG);
    [~, lam] = eig(sqrtG*ERROR.Ae*invsqrtG + invsqrtG*ERROR.Ae'*sqrtG);
    alpha = max(real(diag(lam)))/2;
    if alpha >= 0
        fprintf('alpha is greater than 0, try using a different method.\n');
        return;
    else
        beta = 1;
        sqrtG = sqrtm(G);
    end
    
elseif strcmp(G_method, 'Lyap')
    fprintf('Computing G, beta, alpha via Lyapunov method.\n');
    [~, lam] = eig(full(ERROR.Ae));
    perturb = 0.1;
    success = false;
    while ~success
        try
            alpha = max(real(diag(lam))) + perturb;
            R = lyapchol((ERROR.Ae - alpha*eye(ne))', eye(ne)); % seems to perform better numerically than lyap
            success = true;
        catch
            perturb = perturb*2;
        end
    end
    G = R'*R;
    sqrtG = sqrtm(G);
    invsqrtG = inv(sqrtG);
    beta = 1;

elseif strcmp(G_method, 'LMI')   
    fprintf('Computing G, beta, alpha via LMI constrained optimization.\n');
    [G, alpha, ~] = computeG_LMI_Continuous(ERROR.Ae, opt);
    sqrtG = sqrtm(G);
    beta = 1;

elseif strcmp(G_method, 'GP')
    fprintf('Computing G, beta, alpha via geometric programming.\n');
    [xi, lam] = eig(full(ERROR.Ae));
    invxi = inv(xi);
    alpha = max(real(diag(lam)));
    [G, ~] = computeG_GP({xi, blkdiag(ERROR.Be, ERROR.Ge)}, {invxi, [ERROR.Ez; ERROR.Eu]}, opt);
    sqrtG = sqrtm(G);
    beta = cond(sqrtG*xi);
    
elseif strcmp(G_method, 'BatchGP')
    if all(size(G_0) == [ne, ne]) && isdiag(G_0)
        G_0 = diag(G_0);
    elseif all(size(G_0) == [ne, ne])
        fprintf('G_0 needs to be diagonal to use with BatchGP method.\n');
        return;
    else
        fprintf('G_0 incorrect dimension for BatchGP method, using I to warm start.\n');
        G_0 = ones(ne, 1);
    end
    
    if any(G_0 <= 0)
        fprintf('G_0 needs to be positive definite.\n');
        return;
    end
    
    fprintf('Computing G, beta, alpha via batch geometric programming.\n');
    [xi, lam] = eig(full(ERROR.Ae));
    invxi = inv(xi);
    alpha = max(real(diag(lam)));
    [G] = computeG_BatchGP({xi, blkdiag(ERROR.Be, ERROR.Ge)}, {invxi, [ERROR.Ez; ERROR.Eu]}, ...
                            min(ne, 20), 'random', 10000, 0.001, G_0, opt);
    sqrtG = sqrtm(G);
    
    % This code since I was getting some wierd errors from using cond()
    try
        beta = cond(sqrtG*xi);
    catch
        fprintf('Computing beta1.\n');
        beta1 = matrix2Norm(sqrtG*xi);
        fprintf('Computing beta2.\n');
        beta2 = matrix2Norm(invxi*inv(sqrtG));
        beta = beta1*beta2;
    end
    
else
    if all(size(G_0) == [ne, ne])
        fprintf('Not computing G, using G_0 as specified. Checking if G_0 is PD.\n');
        % Check to make sure G_0 is positive definite
        [~, p] = chol(G_0);
        if p ~= 0
            fprintf('G_0 is not PD.\n');
            return;
        end
        G = G_0;
    else
        fprintf('G_0 incorrect dimension, using G_0 = I instead.\n');
        G = eye(ne);
    end
    
    fprintf('Computing G, beta, alpha with G = G_0.\n');
    sqrtG = sqrtm(G);
    [xi, lam] = eig(full(ERROR.Ae));
    alpha = max(real(diag(lam)));
    beta = cond(sqrtG*xi);
end

% Save the data
if savedata
    G = sparse(G);
    sqrtG = sparse(sqrtG);
    fprintf('Saving data to %s.\n\n', opt.normbound_datapath);
    save(opt.normbound_datapath, 'beta','alpha','G','sqrtG');
end

end

