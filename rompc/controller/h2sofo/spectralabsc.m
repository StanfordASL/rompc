function [f, G] = spectralabsc(A, continuous)

if any(any(isnan(A)|isinf(A)))
    f = inf;
    G = nan(size(A));
    return
end
if nargout < 2
    if continuous
        f = max(real(eig(A)));
    else
        f = max(abs(eig(A)));
    end
else
    [V, Lambda] = eig(A); 
    lambda = diag(Lambda);
    if continuous
        [f, k] = max(real(lambda)); % spectral abscissa
    else
        [f, k] = max(abs(lambda)); % spectral radius
    end
    lam = lambda(k); % corresponding eigenvalue
    v = V(:,k);  % right eigenvector
    I = eye(length(A));
    e = I(k,:);
    warning('off')
    u = e/V;   % relevant row of inverse of V ("left row eigenvector")
    
    % Might possibly have to perturb V to make it nonsingular
    perturb = 1e-16*max(max(abs(V)));
    while isnan(u)
        V = V + perturb*randn(size(V));
        u = e/V;
        perturb = 10*perturb;
    end
    warning('on')
    
    % Compute gradient
    if continuous
        G = u'*v';  % gradient in complex matrix space
    else
        G = (lam/abs(lam))*u'*v';
    end
end