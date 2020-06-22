function [init] = procInitialController(init, pars, options)

% Verify K0 is appropriate dimension
if ~isempty(init) && any(size(init) ~= [pars.m, pars.p])
    error('h2sofo: options.init must have dimension %d by %d', pars.m, pars.p)
end

% Verify that the inputted controller has the proper structure
if ~isempty(init) && isstruct(options.structure) && isfield(options.structure,'Khat') && ~isempty(options.structure.Khat)
    
    % A weird transformation to see what variables are free
    K = options.structure.Khat'*ones(size(options.structure.Khat, 1), 1);
     
    % Check element by element if the fixed parameters are respected in
    % the initial controller
    r = not(logical(reshape(K, size(init))));
    if max(abs(init(r) - options.structure.Kfix(r))) > 0
        error('The initial controller matrix K differs from the user-provided structure');
    end
end

end
