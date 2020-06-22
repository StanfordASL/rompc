function [Acl, Bcl, Ccl, Dcl] = getABCDclosedloop(P, k, structure)

m = size(P.B2, 2);
p = size(P.C2, 1);
[Khat] = k2K(m, p, k, structure);

Acl = P.A + P.B2*Khat*P.C2; % n by n
Bcl = P.B1 + P.B2*Khat*P.D21; % n by m1
Ccl = P.C1 + P.D12*Khat*P.C2; % p1 by n   
Dcl = P.D11 + P.D12*Khat*P.D21; % p1 by m1

end
