function [K] = k2K(m, p, x, structure)

% if K exists it includes the information from Khat
if isfield(structure, 'K')
    if ~isempty(structure.K)
        % unique solution for Khat
        K = structure.K;
    else
        % infinity solution for Khat
        indx = 1:size(structure.V, 2);
        y = x(indx);
        k = structure.V*y + structure.w;
        K = reshape(k, m, p);
    end
elseif isfield(structure, 'Khat')
    % not H2 problem
    if(length(x) ~= size(structure.Khat, 1) )
        fprintf('k2K:length mismatch with structure\n');
    end
    
    K = reshape(structure.Khat'*x, m, p) + structure.Kfix;
else
    %no constraints so just reshape
    indx = 1:m*p;
    if length(x) ~= indx(length(indx))
        fprintf('k2K: length mismatch')
    end
    K = reshape(x(indx), m, p);
end
    
end
