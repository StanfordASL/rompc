function [FOM, ROM, CTRL, Z, U, NOISE, EBOUND, PARAMS] = tubularReactor(regen)
%[FOM, ROM, CTRL, Z, U, NOISE, EBOUND, PARAMS] = tubularReactor(regen)
%
%Loads the system FOM, ROM, controllers, constraints, and disturbances. Defines the state and
%control constraints. Loads error bounds if they have already been
%computed.
%
%Tubular reactor for a chemical process, described by Agudelo et al.
%(2007). The state x = [C_1, ..., C_N, T_1, ..., T_N] where C is the
%concentration at a point along the reactor axis and T is the temperature.
%Both quantities are normalized with respect to a reference temperature.
%The concentration and temperature are in units of C [mol/L] and T [K] when
%dimensionalized.
%
% Inputs:
%   regen: false to just load everything from data files

path = fileparts(mfilename('fullpath'));

% Parameters
P.N = 300; % spatial discretization parameter
P.v = 0.1;
P.L = 1;
P.dz = P.L/P.N;
P.Tnorm = 340;
P.Cnorm = 0.02;
P.k0 = 10^6;
P.E = 11250;
P.R = 1.986;
P.Cin = 0.02;
P.Tin = 340;
P.Gr = 4.25*10^9;
P.Hr = 0.2;

% Compute steady state operating point to define optimal operation
TJ1 = 374.6;
TJ2 = 310.1;
TJ3 = 325.2;
[Cstar, Tstar] = reactorSteadyState(TJ1, TJ2, TJ3, P);

% Define FOM
if regen
    fprintf('Regenerating system data.\n');
       
    % System matrices
    nf = 2*P.N;
    Af = zeros(nf);
    Bf = zeros(nf, 3);
    Bfw = zeros(nf, 2);
    
    % Compute alpha parameters
    alpha_A = P.k0*exp(-P.E./(P.R*Tstar));
    alpha_B = P.k0*P.E*Cstar.*exp(-P.E./(P.R*Tstar))./(P.R*Tstar.^2);
    alpha_C = -P.Gr*exp(-P.E./(P.R*Tstar));
    alpha_D = -P.Gr*P.E*Cstar.*exp(-P.E./(P.R*Tstar))./(P.R*Tstar.^2) + P.Hr;
    
    % Add equations for concentration along reactor axis
    Bfw(1,1) = P.v/P.dz;
    for i = 1:P.N  
        Af(i,i) = -P.v/P.dz - alpha_A(i);
        Af(i,P.N+i) = -alpha_B(i)*P.Tnorm/P.Cnorm;
        if i > 1
            Af(i,i-1) = P.v/P.dz;
        end
    end
    
    % Add equations for temperature along reactor axis
    Bfw(P.N+1,2) = P.v/P.dz;
    for i = P.N+1:nf
        Af(i,i) = -P.v/P.dz - alpha_D(i-P.N);
        Af(i,i-P.N) = -alpha_C(i-P.N)*P.Cnorm/P.Tnorm;
        if i > P.N+1
            Af(i,i-1) = P.v/P.dz;
        end
        
        % Control inputs
        if i - P.N <= 100
            Bf(i,1) = P.Hr;
        elseif i - P.N > 100 && i - P.N <= 200
            Bf(i,2) = P.Hr;
        elseif i - P.N > 200
            Bf(i,3) = P.Hr;
        end
    end
    
    % Output matrix: temperature measurement at certain points
    Cf = zeros(3, nf);
    Cf(1, P.N + 50) = 1;
    Cf(2, P.N + 150) = 1;
    Cf(3, P.N + 250) = 1;

    % Performance matrix: temperature at some points
    Hf = zeros(10, nf);
    for i = 1:10
        Hf(i,P.N + 30*(i-1) + 1) = 1; % N + 1, N + 31, ...
    end
         
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
Tmax = 395;
Tmin = 300;
umax = 395;
umin = 300;
zUB_full = (Tmax - Tstar)/P.Tnorm;
zLB_full = (Tmin - Tstar)/P.Tnorm;
zUB = zUB_full(1:30:end);
zLB = zLB_full(1:30:end);
uUB = [(umax - TJ1)/P.Tnorm, (umax - TJ2)/P.Tnorm, (umax - TJ3)/P.Tnorm];
uLB = [(umin - TJ1)/P.Tnorm, (umin - TJ2)/P.Tnorm, (umin - TJ3)/P.Tnorm];
Z = rectanglePolytope(zUB, zLB);
U = rectanglePolytope(uUB, uLB);

% Create bounds on noise w, v
NOISE.W = 0;
NOISE.V = 0;

% Cost 
FOM.Qf = 1*eye(nf);
FOM.Rf = 100*eye(m);
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

% Parameters that may be useful to output
PARAMS.Hf_full = [zeros(P.N), eye(P.N)];
PARAMS.Z_full = rectanglePolytope(zUB_full, zLB_full);
PARAMS.P = P;
PARAMS.ustar = [TJ1; TJ2; TJ3];
PARAMS.Tstar = Tstar;
PARAMS.Cstar = Cstar;
end

