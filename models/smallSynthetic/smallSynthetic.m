function [FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = smallSynthetic(regen)
%[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = smallSynthetic(regen)
%
%Loads the system FOM, ROM, controllers, constraints, and disturbances. Defines the state and
%control constraints. Loads error bounds if they have already been
%computed.
%
% Inputs:
%   regen: false to just load everything from data files

path = fileparts(mfilename('fullpath'));

% Define FOM
if regen
    fprintf('Regenerating system data.\n');
    
    % System used in Lorenzetti et al. 2019 (ECC)
    Af = [0.28, 0.25, -0.19, -0.22, 0.03, -0.50;
     0.25, -0.47, 0.30, 0.17, -0.11, -0.11;
     -0.19, 0.30, 0.46, 0.09, -0.02, -0.08;
     -0.22, 0.17, 0.09, 0.60, -0.06, 0.14;
     0.03, -0.11, -0.02, -0.06, 0.46, -0.13;
     -0.50, -0.11, -0.08, 0.14, -0.13, -0.23];

    Bf = [1.0159; 0; 0.5988; 1.8641; 0; -1.2155];

    Bfw = eye(6);

    Cf = [1.2920, 0.2361, 0, 0, 0, 0];

    Hf = [1, 0, 0, 0, 0, 0;
          0, 1, 0, 0, 0, 0];
      
    FOM = struct('Af', Af, 'Bf', Bf, 'Bfw', Bfw, 'Cf', Cf, 'Hf', Hf);
      
    % Save data
    fprintf('Saving system data to %s.\n', strcat(path, '/data/FOM.mat'));
    save(strcat(path, '/data/FOM.mat'), 'FOM');
else
    fprintf('Loading system data from %s.\n', strcat(path, '/data/FOM.mat'));
    load(strcat(path, '/data/FOM.mat'));
end

nf = size(FOM.Af, 1);
m = size(FOM.Bf, 2);
p = size(FOM.Cf, 1);
o = size(FOM.Hf,1);
    
% Load ROM
rom_path = strcat(path, '/data/ROM.mat');
if exist(rom_path, 'file')
    fprintf('Loading ROM from %s.\n', rom_path);
    load(rom_path);
else
    ROM = struct();
end

% Load controller
ctrl_path = strcat(path, '/data/CTRL.mat');
if exist(ctrl_path, 'file')
    fprintf('Loading controller from %s.\n', ctrl_path);
    load(ctrl_path);
else
    CTRL = struct();
end

% Constraints
zUB = [50, 50];
zLB = -zUB;
uUB = 20;
uLB = -uUB;
Z = rectanglePolytope(zUB, zLB);
U = rectanglePolytope(uUB, uLB);

% Create bounds on noise w, v
wUB = 0.05*ones(1, nf);
wLB = -0.05*ones(1, nf);
vUB = 0.01*ones(1, p);
vLB = -0.01*ones(1, p);
NOISE.W = rectanglePolytope(wUB, wLB);
NOISE.V = rectanglePolytope(vUB, vLB);

% Cost 
FOM.Qf = 10*eye(nf);
FOM.Rf = 1*eye(m);
if isfield(ROM, 'V')
    ROM.Q = ROM.V'*FOM.Qf*ROM.V; 
    ROM.R = FOM.Rf;
end

% Load error
error_path = strcat(path, '/data/EBOUND.mat');
if exist(error_path, 'file')
    fprintf('Loading error data from %s.\n', error_path);
    load(error_path);
else
    EBOUND = struct();
end

end

