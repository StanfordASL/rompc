clear; clc;
path = fileparts(which('template.m'));

[~, ~, ~, Z, U] = template(false);

% Analyze error bounds WITH noise
load(strcat(path, '/data/EBOUND.mat'));
[Dz, Du, ~, ~] = tightenConstraints(Z, U, EBOUND, 0);
fprintf('z bounds, %% (with noise): ');
for i = 1:length(Dz)
    fprintf('%.1f, ', 100*Dz(i)/Z.b(i));
end
fprintf('\n');
fprintf('u bounds, %% (with noise): ');
for i = 1:length(Du)
    fprintf('%.1f, ', 100*Du(i)/U.b(i));
end
fprintf('\n');

% Analyze error bounds WITHOUT noise
load(strcat(path, '/data/EBOUND_no_noise.mat'));
[Dz, Du, ~, ~] = tightenConstraints(Z, U, EBOUND, 0);
fprintf('z bounds, %% (without noise): ');
for i = 1:length(Dz)
    fprintf('%.1f, ', 100*Dz(i)/Z.b(i));
end
fprintf('\n');
fprintf('u bounds, %% (without noise): ');
for i = 1:length(Du)
    fprintf('%.1f, ', 100*Du(i)/U.b(i));
end
fprintf('\n');