clear; clc; close all;

A = [1, -1; 1, 1];
B = ones(2,1);
C = [1, 0];
H = C;
Bw = ones(2,1);
P=struct('A',A,'B1',Bw,'B2',B,...
'C1',H,'C2',C,'D11',zeros(size(H,1), size(Bw,2)),...
'D12',zeros(size(H,1), size(B,2)),'D21',zeros(size(C,1), size(Bw,2)));

options.fast = 1;
options.prtlevel=2;
options.discrete = true;
options.init = [0];
[K, f] = h2sofo(P, 'h2', options);

% Verify optimum
N = 1000;
Ks = linspace(-2,2,N);
fs = zeros(size(Ks));
for i = 1:N
    Acl = P.A + P.B2*Ks(i)*P.C2; % n by n
    Bcl = P.B1 + P.B2*Ks(i)*P.D21; % n by m1
    Ccl = P.C1 + P.D12*Ks(i)*P.C2; % p1 by n   

    % Check stability
    [V,D] = eig(Acl);
    flam = max(abs(diag(D)));
    if flam < 1
        % Compute the H2 norm
        [X] = dlyap(Acl', Ccl'*Ccl);
        fs(i) = sqrt(trace(Bcl'*X*Bcl));
    else
        fs(i) = nan;
    end
end
figure; hold on;
plot(Ks,fs);
plot(K, f, 'r*');




