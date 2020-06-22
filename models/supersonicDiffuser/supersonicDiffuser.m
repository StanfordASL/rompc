function [FOM, ROM, CTRL, Z, U, NOISE, EBOUND, params] = supersonicDiffuser()
%[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = supersonicDiffuser()
%
%Loads the system FOM, ROM, controllers, constraints, and disturbances. Defines the state and
%control constraints. Loads error bounds if they have already been
%computed.
%
% Inputs:
%   regen: false to just load everything from data files

path = fileparts(mfilename('fullpath'));

fprintf('Regenerating system data.\n');
    
% Load data
load(strcat(path, '/definition/inlet.mat'));
Af = Problem.A;
Bf = Problem.aux.B;
Ef = Problem.aux.E;
Cf = Problem.aux.C;
Hf = Cf;

% Discretize using backward Euler
dt = .025;
D = Ef + dt*Af;
invD = inv(full(D));
Af = invD*Ef;
Bf = dt*invD*Bf;
   
% Separate out control and disturbance
Bfw = Bf(:,2);
Bf = Bf(:,1);

% Now convert to a system where the disturbance is a state
% such that the state is [xf; d] and the disturbance is the change
% in d
Af = [Af, Bfw; zeros(1, size(Problem.A,1)), 0.95];
Bf = [Bf; 0];
Bfw = [zeros(size(Problem.A,1),1); 1];
Cf = [Cf, 0];
Hf = [Hf, 0; zeros(1,size(Problem.A,1)), 1]; % z = [z; d]
      
FOM = struct('Af', Af, 'Bf', Bf, 'Bfw', Bfw, 'Cf', Cf, 'Hf', Hf);

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

% Parameters for w = M*exp(-alpha*(t-t0)^2)
f0 = 0.1; % since dt = 1/4 f0 and dt = 0.25
params.dt = dt;
params.M = 0.01*1.225; % equal to 1% sea level density
params.alpha = 2*f0^2; % same choice as Gravedahl et al.

% Constraints
% yss is about 1.36 (throat Mach number)
ymax = 0.44; % keep throat Mach number below 1.8
ymin = -0.26; % keep throat Mach number greater than 1.1
bmax = 0.04; % max bleed 
bmin = -0.01; % keeps total bleed positive (steady state + pert)
zUB = [ymax, params.M];
zLB = [ymin, -params.M];
uUB = [bmax];
uLB = [bmin];
Z = rectanglePolytope(zUB, zLB);
U = rectanglePolytope(uUB, uLB);

% Create bounds on noise w, v
wUB = [0.00065]; % rate of change of the density disturbance
wLB = -wUB;
vUB = [.0001];
vLB = -vUB;
NOISE.W = rectanglePolytope(wUB, wLB);
NOISE.V = rectanglePolytope(vUB, vLB);

% Cost 
FOM.Qf = blkdiag(speye(nf-1), [0]);
FOM.Rf = 0.05*eye(m);
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

