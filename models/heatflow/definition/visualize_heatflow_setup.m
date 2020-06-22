clear; clc; close all;
path = fileparts(which('visualize_heatflow_setup.m'));
load(strcat(path, '/hf2d9.mat'));
[Af, Bf, Cf, Hf] = heatflowMods(Af);


% Plot results
figure(); hold on;
labels = {};
for i = 1:size(Bf,2)
    [T, posx, posy] = createTempMesh(Bf(:,i));
    idx = find(T>0);
    plot(posx(idx), posy(idx), 'b.');
    labels{i} = sprintf('Control %d', i);
end

% Performance locations
for i = 1:size(Hf,1)
    [T, posx, posy] = createTempMesh(Hf(i,:)');
    idx = find(T>0);
    plot(posx(idx), posy(idx), 'r.');
    labels{size(Bf,2) + i} = sprintf('Performance %d', i);
end
legend(labels);
xlim([0 60]);
ylim([0 60]);
axis square