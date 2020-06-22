function [x,N,x0] = solveSystem(A,b)
%   solves the linear sistem Ax = b
%   returns x = nan if no solution, x = [] for infinity solution or x if
%   unique solution
%   for the infinity solution x = x0 + N*y where y is arbitrary

[~,n] = size(A);
[m,p] = size(b);

[U,S,V] = svd(A);   %   A = U*S*V'

%solve the system S*y = c with y = V'*x and c = U'*b

c = U'*b;
r = rank(S,1e-10*S(1,1));

idx = [];
if r < m                    %then we should check the norm in c
    idx = r+1:m;
end

if norm(c(idx)) > 1e-10*S(1,1)        %no solution - threshold chosen randomly
    x = nan;
    N = nan;
    x0 = nan;
    nf = nan;
    return;
else
    nf = n-r;
    y = zeros(r,p);
    % unique solution or infinity solution
    for i = 1:r
        y(i) = c(i)/S(i,i);
    end
        
    if (nf == 0)            %unique solution
        x = V*y;
        x0 = [];
        N = [];
    else
        %infinity solution
        x0 = V*[y; zeros(nf,1)];
        x = [];
        N = V*[zeros(r,nf); eye(nf)];
    end
end
