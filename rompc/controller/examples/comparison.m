clear; clc; close all;

%% Discrete time Random small system
opts = {};
data = {};
[FOM, ROM, ~, ~, ~, ~, ~] = randomSmallSystem(false);
opts.discrete = true;
opts.with_noise = false;

% Compute LQR controller
opts.method = 'ric';
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.ric_time = toc(t0);
CTRL_Ric.K = K;
CTRL_Ric.L = L;
[data.ric] = computeErrorH2Norm(FOM, ROM, CTRL_Ric, opts);

% Compute H2SOFO controller with warmstart
opts.method = 'h2sofo';
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.hifoo_time = toc(t0);
CTRL_H2SOFO.K = K;
CTRL_H2SOFO.L = L;
[data.hifoo] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO, opts);

% Compute H2SOFO controller without warmstart
opts.method = 'h2sofo';
opts.warmstart = true;
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.hifoo_ric_time = toc(t0);
CTRL_H2SOFO.K = K;
CTRL_H2SOFO.L = L;
[data.hifoo_ric] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO, opts);

fprintf('Riccati: %.4f, H2SOFO: %.4f, H2SOFO+Riccati: %.4f\n', data.ric, data.hifoo, data.hifoo_ric);
fpath = strcat(fileparts(which('h2optController.m')), '/examples');
save(strcat(fpath, '/randomSmall_compare.mat'), 'data');

%% Continuous time tubular reactor
opts = {};
data = {};
[FOM, ROM, ~, ~, ~, ~, ~, ~] = tubularReactor(false);
opts.continuous = true;
opts.with_noise = false;

% Compute LQR controller
opts.method = 'ric';
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.ric_time = toc(t0);
CTRL_Ric.K = K;
CTRL_Ric.L = L;
[data.ric] = computeErrorH2Norm(FOM, ROM, CTRL_Ric, opts);

% Compute H2SOFO controller with warmstart
opts.method = 'h2sofo';
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.hifoo_time = toc(t0);
CTRL_H2SOFO.K = K;
CTRL_H2SOFO.L = L;
[data.hifoo] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO, opts);

% Compute H2SOFO controller without warmstart
opts.method = 'h2sofo';
opts.warmstart = true;
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.hifoo_ric_time = toc(t0);
CTRL_H2SOFO.K = K;
CTRL_H2SOFO.L = L;
[data.hifoo_ric] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO, opts);

fprintf('Riccati: %.4f, H2SOFO: %.4f, H2SOFO+Riccati: %.4f\n', data.ric, data.hifoo, data.hifoo_ric);
fpath = strcat(fileparts(which('h2optController.m')), '/examples');
save(strcat(fpath, '/tubularReactor_compare.mat'), 'data');


%% Discrete time heat flow
opts = {};
data = {};
[FOM, ROM, ~, ~, ~] = distHeatflowDiscrete();
opts.discrete = true;
opts.with_noise = false;

% Compute LQR controller
opts.method = 'ric';
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.ric_time = toc(t0);
CTRL_Ric.K = K;
CTRL_Ric.L = L;
[data.ric] = computeErrorH2Norm(FOM, ROM, CTRL_Ric, opts);

% Compute H2SOFO controller with warmstart
opts.method = 'h2sofo';
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.hifoo_time = toc(t0);
CTRL_H2SOFO.K = K;
CTRL_H2SOFO.L = L;
[data.hifoo] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO, opts);

% Compute H2SOFO controller without warmstart
opts.method = 'h2sofo';
opts.warmstart = true;
t0 = tic;
[K, L] = h2optController(FOM, ROM, opts);
data.hifoo_ric_time = toc(t0);
CTRL_H2SOFO.K = K;
CTRL_H2SOFO.L = L;
[data.hifoo_ric] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO, opts);

fprintf('Riccati: %.4f, H2SOFO: %.4f, H2SOFO+Riccati: %.4f\n', data.ric, data.hifoo, data.hifoo_ric);
fpath = strcat(fileparts(which('h2optController.m')), '/examples');
save(strcat(fpath, '/distHeatflow_compare.mat'), 'data');

