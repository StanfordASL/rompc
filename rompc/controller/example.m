clear; clc; close all;

[FOM, ROM] = smallSynthetic(false);

opts.discrete = true;
opts.with_noise = false;

% Compute controller via reduced order Riccati method
opts.method = 'ric';
[CTRL_RIC.K, CTRL_RIC.L] = computeControllerGains(FOM, ROM, opts);
[n_ric] = computeErrorH2Norm(FOM, ROM, CTRL_RIC, opts);

% Compute controller using gradient based optimization, with warmstart from
% reduced order Riccati method
opts.method = 'h2sofo';
opts.warmstart = true;
[CTRL_H2SOFO_RIC.K, CTRL_H2SOFO_RIC.L] = computeControllerGains(FOM, ROM, opts);
[n_H2SOFO_RIC] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO_RIC, opts);

% Compute controller using gradient based optimization, with warmstart from
% reduced order Riccati method
opts.method = 'h2sofo';
opts.warmstart = false;
[CTRL_H2SOFO.K, CTRL_H2SOFO.L] = computeControllerGains(FOM, ROM, opts);
[n_H2SOFO] = computeErrorH2Norm(FOM, ROM, CTRL_H2SOFO, opts);

fprintf('Riccati: %.4f, H2SOFO+Riccati: %.4f, H2SOFO: %.4f\n', n_ric, n_H2SOFO_RIC, n_H2SOFO);

