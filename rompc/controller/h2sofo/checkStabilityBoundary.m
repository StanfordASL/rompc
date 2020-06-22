function checkStabilityBoundary(ksol, pars, prtlevel)
% Verify that the eigenvalues of the closed loop A
% matrix are remaining within the stability boundary
% Let x+yi be the eigenvalues of Acl.  We will
% ensure that x/(1+x^2+y^2) <= 1e-6

[Acl, ~, ~, ~] = getABCDclosedloop(pars.plantinfo, ksol, pars.structure);

evals = eig(Acl);

if pars.continuous
    boundary = real(evals)./(1+real(evals).^2+imag(evals).^2);
    if nnz(boundary >= -1e-6) > 0
        if prtlevel > 0
            fprintf('h2sofo: An eigenvalue of the closed-loop system is approaching the stability boundary\n');
        end
    end
else
    boundary = abs(evals);
    if nnz(boundary >= 1 - 1e-6) > 0
        if prtlevel > 0
            fprintf('h2sofo: An eigenvalue of the closed-loop system is approaching the stability boundary\n');
        end
    end
end
 
end