function [Kstr, Kfix, Kno] = getStructure(Kstruct)
%
% Inputs:
%   Kstruct: the matrix defining the structure, with nan for opt variables
%
% Returns:
%   Kstr: the new structure matrix obtained by selecting from the I matrix
%   Kfix: matrix of the same size as Kstruct containing fixed params and 0s
%   Kno: specifies which elements should not be optimized 


[m, p] = size(Kstruct);
npars = m*p;

I = eye(m*p);
vect = Kstruct(:);
Kstr = []; 
Kno = []; 
Kfix = vect;
for i = 1:npars
    if isnan(vect(i)) % parameter should be optimized
        Kstr = [Kstr; I(i,:)];
        Kfix(i) = 0;
    else
        Kno = [Kno; I(i,:)];
    end
end
Kfix = reshape(Kfix, m, p);

end
