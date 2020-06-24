function [M, gamma, G, sqrtG] = computeScalarDecayBoundDiscrete(ERROR, G_method, G_0, opt)
%[M, gamma, G, sqrtG] = computeScalarDecayBoundDiscrete(ERROR, G_method, G_0, opt)
%
%Computes a bound ||Ae^k||_G \leq M*gamma^k for discrete time problems
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
%   M: constant in norm bound
%   gamma: exponential term constant in norm bound
%   G: PD matrix defining weighted norm x^TGx = ||x||_G^2
%   sqrtG: matrix square root of G

ne = size(ERROR.Ae, 1);

% First check if data can be loaded, or if it should be saved
if isfield(opt, 'normbound_datapath')
    if exist(opt.normbound_datapath, 'file')
        fprintf('******************************\n');
        fprintf('Loading norm transient bound gamma, M from %s.\nDelete this file to recompute.\n', opt.normbound_datapath);
        fprintf('******************************\n\n');
        load(opt.normbound_datapath, 'M', 'gamma', 'G', 'sqrtG');
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

if strcmp(G_method, 'Lyap')
    fprintf('Computing G, M, gamma via discrete Lyapunov method.\n');
    [~, lam] = eig(full(ERROR.Ae));
    gamma = max(abs(diag(lam))) + .0001;
    G = dlyap(ERROR.Ae'/gamma, eye(ne)/gamma^2);
    sqrtG = sqrtm(G);
    gamma = norm(sqrtG*ERROR.Ae*inv(sqrtG));
    M = 1;

elseif strcmp(G_method, 'LMI')
    fprintf('Computing G, M, gamma via LMI constrained optimization.\n');
    [G, ~] = computeG_LMI_Discrete(ERROR.Ae, opt);
    sqrtG = sqrtm(G);
    gamma = norm(sqrtG*ERROR.Ae*inv(sqrtG));
    M = 1;

elseif strcmp(G_method, 'GP')
    fprintf('Computing G, M, gamma via geometric programming.\n');
    [xi, lam] = eig(full(ERROR.Ae));
    invxi = inv(xi);
    gamma = max(abs(diag(lam)));
    [G, ~] = computeG_GP({xi, blkdiag(ERROR.Be, ERROR.Ge)}, {invxi, [ERROR.Ez; ERROR.Eu]}, opt);
    sqrtG = sqrtm(G);
    M = cond(sqrtG*xi);
    
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
    
    fprintf('Computing G, M, gamma via batch geometric programming.\n');
    [xi, lam] = eig(full(ERROR.Ae));
    invxi = inv(xi);
    gamma = max(abs(diag(lam)));
    [G] = computeG_BatchGP({xi, blkdiag(ERROR.Be, ERROR.Ge)}, {invxi, [ERROR.Ez; ERROR.Eu]}, ...
                            min(ne, 20), 'random', 10000, 0.001, G_0, opt);
    sqrtG = sqrtm(G);
    
    % This code since I was getting some wierd errors from using cond()
    try
        M = cond(sqrtG*xi);
    catch
        fprintf('Computing M1.\n');
        M1 = matrix2Norm(sqrtG*xi);
        fprintf('Computing M2.\n');
        M2 = matrix2Norm(invxi*inv(sqrtG));
        M = M1*M2;
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
    
    fprintf('Computing G, M, gamma with G = G_0.\n');
    sqrtG = sqrtm(G);
    [xi, lam] = eig(full(ERROR.Ae));
    gamma = max(abs(diag(lam)));
    M = cond(sqrtG*xi);
end

% Save the data
if savedata
    G = sparse(G);
    sqrtG = sparse(sqrtG);
    fprintf('Saving data to %s.\n\n', opt.normbound_datapath);
    save(opt.normbound_datapath, 'M','gamma','G','sqrtG');
end

end

