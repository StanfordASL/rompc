function g = chainK(P, GA, structure)
% Construct the gradient with respect to the parameterization specified
% in function getABCDbig, via chain rule.
% The complex gradients of the objective function of Abig, Bbig, Cbig, Dbig
% are given by GA, GB, GC, GD respectively. In principle, the objective
% function could be any function, but we are interested in two cases: the
% H-infinity norm (see hinfty.m) and the spectral abscissa (see specabsc.m;
% in this case GB, GC and GD are all trivially zero).

% Derivation of chain rule is from first principles: 
% f(Q + M(X+DeltaX) N) is approx f(Q  + MXN) + <del f, M DeltaX N> and the 
% latter is real tr G' M DeltaX N = real tr N G' M DeltaX = <M' G N', DeltaX>,
% where G is the gradient of f in matrix space.  Thus the gradient wrt X
% of f(Q + MXN) is M'GN'; equivalently, the gradient wrt vec(X) is vec(M'GN').


% Note: the reason we take the real part is not because of the real inner
% product on complex matrix space.  It's because the parameter space is always
% real, even though the matrix functions may be complex.  Usually the matrix
% functions are also real, but of course the eigenvalues are generally complex.

indx = 1:size(P.A,1);
delKhat = real(P.B2'*GA(indx,indx)*P.C2'); % [B2' 0]*GA*[C2'; 0]                                        
g = K2k(delKhat, structure);

end
    
    
    
    
    
    
    
    