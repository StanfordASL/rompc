function [Af, Bf, Cf, Hf] = heatflowMods(Af)

nf = size(Af,1);
N = sqrt(nf);

% Create new control points
patchsize = [5, 5, 5, 5, 5];
upleft_x = [28, 1, 28, 55, 28]; % with 1 being far left of screen and 59 - patchsize far right
upleft_y = [1, 28, 28, 28, 55]; % with 1 being top
m = length(upleft_x);
Bf = zeros(nf,m);
for k = 1:m
    Bf(:,k) = sparse(size(Af,1),1);
    for i = 1:patchsize(k)
        for j = 1:patchsize(k)
            Bf((upleft_x(k) - 2 + i)*N + upleft_y(k) - 1 + j, k) = 1;
        end
    end
end

% Create measurement points
Cf = Bf';
for i = 1:size(Cf,1)
    Cf(i,:) = Cf(i,:)/sum(Cf(i,:));
end

% Create performance points
patchsize = [3, 3, 3, 3];
upleft_x = [5, 44, 40, 10]; % with 1 being far left of screen and 59 - patchsize far right
upleft_y = [5, 14, 40, 40]; % with 1 being top
o = length(upleft_x);
Hf = zeros(nf,o);
for k = 1:o
    Hf(:,k) = sparse(size(Af,1),1);
    for i = 1:patchsize(k)
        for j = 1:patchsize(k)
            Hf((upleft_x(k) - 2 + i)*N + upleft_y(k) - 1 + j, k) = 1;
        end
    end
end
Hf = Hf';
for i = 1:size(Hf,1)
    Hf(i,:) = Hf(i,:)/sum(Hf(i,:)); % average temperature
end

% Update system to fix some stability issues
Af = Af - 0.1804*speye(nf);
end

