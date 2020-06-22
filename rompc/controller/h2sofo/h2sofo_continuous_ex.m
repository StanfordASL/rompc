clear; clc;

% F16 example from Millstone 2006
A = [-0.322, 0.064, 0.0364, -0.9917, 0.0003, 0.0008, 0, 0;
    0, 0, 1, 0.0037, 0, 0, 0, 0;
    -30.6492, 0, -3.6784, 0.6646, -0.7333, 0.1315, 0, 0;
    8.5395, 0, -0.0254, -0.4764, -0.0319, -0.062, 0, 0;
    0, 0, 0, 0, -20.2, 0, 0, 0;
    0, 0, 0, 0, 0, -20.2, 0, 0;
    0, 0, 0, 57.2958, 0, 0, -1, 0;
    0, -1, 0, 0, 0, 0, 0, 0];
B2 = zeros(8,2); B2(5,1) = 20.2; B2(6,2) = 20.2;
C2 = [0, 0, 0, 0, 0, 0, 0, 1;
      0, 0, 0, -57.2958, 0, 0, 1, 0;
      0, 0, 1, 0, 0, 0, 0, 0;
      0, -1, 0, 0, 0, 0, 0, 0];
B1 = zeros(8,2); B1(8,1) = 1;
C1 = [0, 1, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 57.2958, 0, 0, -1, 0];
D12 = zeros(2,2); 
D11 = zeros(2,2); 
D21 = [0, 0;
       0, 1;
       0, 0;
       1, 0];

P=struct('A',A,'B1',B1,'B2',B2,...
'C1',C1,'C2',C2,'D11',D11,...
'D12',D12,'D21',D21);
options.continuous = true;
options.prtlevel = 2;
options.init = zeros(2,4);
options.structure = [nan, 0, nan, nan; 
                      0, nan, 0, 0];
[K, f] = h2sofo(P, 'h2', options);

Acl = P.A + P.B2*K*P.C2; % n by n
Bcl = P.B1 + P.B2*K*P.D21; % n by m1
Ccl = P.C1 + P.D12*K*P.C2; % p1 by n   

% Compute the H2 norm
[X] = lyap(Acl', Ccl'*Ccl);
fH2 = sqrt(trace(Bcl'*X*Bcl));

[V,D] = eig(Acl);
flam = max(real(diag(D)));

