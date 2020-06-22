function [FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = distillationColumn(regen)
%[FOM, ROM, CTRL, Z, U, NOISE, EBOUND] = distillationColumn(regen)
%
%Loads the system FOM, ROM, controllers, constraints, and disturbances. Defines the state and
%control constraints. Loads error data if it exists.
%
%From Skogestad, S. and Morari, M. "Understanding the dynamic behavior of
%distillation columns", 1988.
%
% Inputs:
%   regen: false to just load everything from data files

path = fileparts(mfilename('fullpath'));

% Define FOM
if regen
    fprintf('Regenerating system data.\n');
    
    load(strcat(path, '/definition/cola_init.mat'));
    Ls=2.706; Vs=3.206; Ds=0.5; Bs=0.5; Fs=1.0; zFs=0.5;
    [Af, Bf, Hf, ~]=cola_linearize('cola4_lin',Xinit',[Ls Vs Ds Bs Fs zFs]);
    Cf = zeros(10, size(Af,1));
    Cf(1,1) = 1;
    Cf(2,11) = 1;
    Cf(3,21) = 1;
    Cf(4,31) = 1;
    Cf(5,41) = 1;
    Cf(6,41+1) = 1;
    Cf(7,41+11) = 1;
    Cf(8,41+21) = 1;
    Cf(9,41+31) = 1;
    Cf(10,41+41) = 1;

    % Separate inputs and disturbances
    Bfw = Bf(:,5:6);
    Bf = Bf(:,1:4);
    
    % Add controller dynamics ud = delta_u
    Af = [Af, Bf;
          zeros(4,86)];
    Bf = [zeros(82,4); eye(4)];
    Bfw = [Bfw; zeros(4,2)];
    Cf = [Cf, zeros(10,4)];
    Hf = [Hf, zeros(4,4);
          zeros(4,82), eye(4)];
    
    % Discretize
    dt = 1; % minutes
    [Af, Bf, Bfw] = zoh(dt, Af, Bf, Bfw);
    
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
zc = [0.01, 0.01, 1, 1]; % constraints on real z = [y_D, x_B, M_D, M_B]
uc = [1, 1, 1, 1]; % constraints on real u = [L, V, D, B]
duc = [0.2, 0.2, 0.2, 0.2]; % constraints on change in u
zUB = [zc, uc];
zLB = -zUB;
uUB = duc;
uLB = -uUB;
Z = rectanglePolytope(zUB, zLB);
U = rectanglePolytope(uUB, uLB);

% Create bounds on noise w, v
wUB = [0.01, 0.005]; % this is 1% of nominal
wLB = -wUB;
vUB = 0.0001*ones(1,p);
vLB = -vUB;
NOISE.W = rectanglePolytope(wUB, wLB);
NOISE.V = rectanglePolytope(vUB, vLB);

% Cost
Qz = eye(4);
Qu = 0.01*eye(4);
Qdu = 0.01*eye(4);
FOM.Qf = FOM.Hf'*blkdiag(Qz, Qu)*FOM.Hf;
FOM.Rf = Qdu;
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

