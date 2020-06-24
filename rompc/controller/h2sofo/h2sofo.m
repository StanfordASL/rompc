function [K, f, loc, options] = h2sofo(P, obj, options)
%[K, f, loc, options] = h2sofo(P, obj, options)
%
% Solve H2 optimal static output feedback controller synthesis problems
% for continuous or discrete time systems. Optimization is performed by
% using a gradient based optimization algorithm HANSO. This code is a modification
% to the open-source Matlab package HIFOO. For details on HIFOO see
%
% D. Arzelier, G. Deaconu, S. Gumussoy and D. Henrion. H2 for HIFOO
% International Conference on Control and Optimization with Industrial Applications, 
% Bilkent University, Ankara, Turkey, August 2011
%
% and
%
% S. Gumussoy, D. Henrion, M. Millstone and M.L. Overton.
% Multiobjective Robust Control with HIFOO 2.0
% Proceedings of the IFAC Symposium on Robust Control Design, Haifa, 2009
%
%
% It is assumed the dynamics are given as: (or discrete time)
% xd = Ax + B1w + B2u
% z = C1x + D11w + D12u
% y = C2x + D21w
%
% The code searches for a controller u = Ky (where K may have structure)
%
% Inputs:
%   P: struct with problem data A,B1,B2,C1,C2,D11,D12,D21
%   obj: 'h2' for H2 and 's' for spectral abscissa (radius)
%   options: options structure with
%       - continuous (discrete): required, specifies if problem is discrete
%                                or continuous time
%       - init: initial values for controller gains
%       - prtlevel: 0, 1, ... for more or less printed output
%
% Outputs:

if ~isfield(options, 'init')
    options.init = [];
end

if ~isfield(options, 'prtlevel')
    options.prtlevel = 1;
end

if ~isfield(options, 'continuous') && ~isfield(options, 'discrete')
    error('Need to specify continuous or discrete system.\n');
elseif isfield(options, 'continuous')
    pars.continuous = options.continuous;
elseif isfield(options, 'discrete')
    pars.continuous = ~options.discrete;
end

pars.plantinfo = P;
pars.plantinfo.objective = obj;
pars.m = size(P.B2, 2);
pars.p = size(P.C2, 1);

% Set the structure for the controller
[options.structure] = procStructure(pars, options);

% Check that the initial controller has correct structure
procInitialController(options.init, pars, options);

% Process other options
[pars, options] = procOptions(pars, options);

% Solve the problem
[K, f, loc, options] = h2sofomain(pars, options.init, obj, options);
end
