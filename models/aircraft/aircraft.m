function [FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = aircraft()
%[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = aircraft()
%
%Loads the system FOM, ROM, controllers, constraints, and disturbances. Defines the state and
%control constraints. Loads error bounds if they have already been
%computed.
%
%Note the FOM is not included here because of its large size. The model is
%from:
%
%A. McClellan, J. Lorenzetti, M. Pavone, and C. Farhat, 
%Projection-based Model Order Reduction for Flight Dynamics and Model
%Predictive Control, AIAA Scitech Forum (2020)

path = fileparts(mfilename('fullpath'));


Af = [];
Bf = [];
Cf = [];
Bfw = [];
Hf = [];

FOM = struct('Af', Af, 'Bf', Bf, 'Bfw', Bfw, 'Cf', Cf, 'Hf', Hf);

nf = 998936;
m = 2; % thrust and moment
p = 6; % rigid body states
o = 6; % rigid body states
    
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
zUB = [20, 20, deg2rad(10), 10, 5, deg2rad(10)];
zLB = -zUB;
uUB = [100, 50];
uLB = -uUB;
Z = rectanglePolytope(zUB, zLB);
U = rectanglePolytope(uUB, uLB);

% Create bounds on noise w, v
NOISE.W = 0;
NOISE.V = 0;

% Cost 
FOM.Qf = sparse(1:nf,1:nf,[10, 10, 1000, 100, 100, 100, .001*ones(1, nf-6)],nf,nf);
FOM.Rf = eye(m);
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

