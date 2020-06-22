function x = K2x(K, structure)

if isstruct(structure) == 0
    x = K(:);
else
    if (isfield(structure,'K'))
        if ~isempty(structure.K) % unique solution so there is nothing to include in the optimization set
            x = [];
        else
            % infinity solution for Khat, so y = w + Vx for arbitrary x
            y = K(:);
            x = structure.V'*(y - structure.w); % choose just the optimisation vars for H2 constraint 
        end
    elseif isfield(structure,'Khat')
        % not H2 but we have Khat structure constraint
        x = structure.Khat*K(:);
    else
        %no constraints on Dhat
        x = K(:);
    end     
end

end