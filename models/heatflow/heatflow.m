function [FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = heatflow(regen)
%[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = heatflow(regen)
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
    
    % Load the continuous time heatflow model (from COMPleib)
    load(strcat(path, '/definition/hf2d9.mat'), 'Af', 'Bf', 'Bfw', 'Cf', 'Hf');
    
    % Make some modifications to make the problem more interesting
    [Af, Bf, Cf, Hf] = heatflowMods(Af);
    
    % Discretize the system
    fprintf('Converting FOM to discrete time.\n');
    dt = 0.05;
    [Af, Bf, Bfw] = zoh(dt, Af, Bf, Bfw);
      
    FOM = struct('Af', Af, 'Bf', Bf, 'Bfw', Bfw, 'Cf', Cf, 'Hf', Hf, 'dt', dt);
      
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
zUB = 2*ones(1, o);
zLB = -2*ones(1, o);
uUB = 100*ones(1, m);
uLB = -100*ones(1, m);
Z = rectanglePolytope(zUB, zLB);
U = rectanglePolytope(uUB, uLB);

% No noise for this example
NOISE.W = 0;
NOISE.V = 0;

% Cost 
FOM.Qf = sparse(diag(ones(nf,1))); % state cost
FOM.Rf = 0.01*eye(m); % control cost
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

