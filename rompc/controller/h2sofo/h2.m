function [f, g] = h2(P, k, structure, continuous)

% Get the closed loop system with controller parameters k
[Acl, Bcl, Ccl, ~] = getABCDclosedloop(P, k, structure);

% verify that the feedback system is stable
vp = eig(Acl);
if continuous
    absc = max(real(vp));
    if (absc >= 0)
        f = inf;
        if nargout > 1
            g = nan*ones(size(structure.V, 2),1);
        end
        return;
    end
    
    % Compute the H2 norm
    [X] = lyap(Acl', Ccl'*Ccl);

    f = sqrt(trace(Bcl'*X*Bcl));

    if nargout > 1
        % Compute the gradient
        [Y] = lyap(Acl, Bcl*Bcl');

        G = 2*((P.B2'*X + P.D12'*Ccl)*Y*P.C2' + P.B2'*X*Bcl*P.D21');
        G = G/(2*f);

        % We have a reduced set of optimization variables imposed by the structure
        % since for H2 optimization you must have structure.V
        g = structure.V'*G(:);

    end
else
    rad = max(abs(vp));
    if (rad >= 1)
        f = inf;
        if nargout > 1
            g = nan*ones(size(structure.V, 2),1);
        end
        return;
    end
    
    % Compute the H2 norm
    [X] = dlyap(Acl', Ccl'*Ccl);

    f = sqrt(trace(Bcl'*X*Bcl));

    if nargout > 1
        % Compute the gradient
        [Y] = dlyap(Acl, Bcl*Bcl');

        G = 2*((P.B2'*X*Acl + P.D12'*Ccl)*Y*P.C2' + P.B2'*X*Bcl*P.D21');
        G = G/(2*f);

        % We have a reduced set of optimization variables imposed by the structure
        % since for H2 optimization you must have structure.V
        g = structure.V'*G(:);

    end
end

end
    

