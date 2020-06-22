function [structure] = procStructure(pars, options)

if ~isfield(options,'structure')
    structure = [];
else      
    % Check structure of K
    [m, p] = size(options.structure);
    if  any([m, p] ~= [pars.m, pars.p])
        error('Error in structure: K has incorrect dimensions');
    end

    [structure.Khat, structure.Kfix, Kno] = getStructure(options.structure);
end
            
% Construct system to check H2 feasibility (see Arzelier, Georgia, et al. 2010
A = kron(pars.plantinfo.D21', pars.plantinfo.D12);
b = -pars.plantinfo.D11(:);

% If structural constraint are enforced, need to account for those as well
if (isstruct(structure) && ~isempty(Kno) )
    A = [A; Kno];
    b = [b; Kno*structure.Kfix(:)];
end

[x, V, w] = solveSystem(A,b);

% K is the controller gain matrix
if (isnan(x)) % no solution
    K = nan;
    error('Direct feedthrough exists, H2 norm is infinite for all controller gains K');
elseif ~isempty(x) % unique solution
    K = reshape(x, m, p);
    error('A unique K exists to ensure no direct feedthrough, no optimization possible');
else % infinite solutions, parametrized through V and w
    K = [];
end

% The controller is parameterized as x = w + V*y and K = reshape(x, m, p)
% over the new decision variables y (which may be fewer than before, which
% is the whole point)
structure.K = K;
structure.V = V;
structure.w = w;

end
