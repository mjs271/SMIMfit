%===============================================================================
%===============================================================================
%
% This function specifies the options to be provided to the optimization
% function in SMIMfit().
% NOTE: currently the optimization function is hard-coded to be lsqnonlin, but
% this could be generalized, if desired.
%
%===============================================================================
%
% Inputs:
%   mymeas := measurement (conservative or reactive)
%
% Outputs:
%   fopts           := optimization options data type (see Matlab documentation
%                      for optimoptions() or lsqnonlin() for further information)
%   iterFit         := number of iterations of fitting to be conducted
%   iterRandRestart := number of iterations for a local random restart of
%                      the optimization problem to help avoid local minima
%   varargout       := optional output (only when called during the reactive
%                      fit), holding:
%       - fitTransParams: a logical flag indicating whether to also fit
%                         transport parameters when conducting the reactive fit
%                         (i.e., "false" implies only reactive parameters will
%                          be passed to the objective function for fitting, and
%                          the fit will run faster)
%
%===============================================================================
%===============================================================================

function [fopts, iterFit, iterRandRestart, varargout] = getOptInputs(mymeas)

p = inputParser;

% error checking function
checkMeas = @(x) any(validatestring(x, {'c', 'r'}));

% measurement, must be 'c' (conservative) or 'r' (reactive)
addRequired(p, 'mymeas', checkMeas);

parse(p, mymeas);

% inputs for the conservative case
if strcmp(mymeas, 'c')
    fitFunc = @lsqnonlin;
    dispOpt = 'iter'; % displays output after each iteration
    funcTol = 5e-8;
    stepTol = 1e-7;
    optTol  = 1e-7;

%     number of iterations for optimization fit
    iterFit = 3;
%     number of iterations for fit with random restart
    iterRandRestart = 0;

    fopts = optimoptions(fitFunc, 'Display', dispOpt, 'FunctionTolerance',...
                         funcTol, 'StepTolerance', stepTol,...
                         'OptimalityTolerance', optTol);

% inputs for the reactive case
elseif strcmp(mymeas, 'r')
    fitFunc = @lsqnonlin;
    dispOpt = 'iter'; % displays output after each iteration
    funcTol = 5e-9;
    stepTol = 1e-9;
    optTol  = 1e-8;

%     logical flag for whether to fit transport parameters, along with
%     reactive parameters (i.e., "false" implies only reactive parameters
%     will be passed to the objective function for fitting)
    fitTransParams = false;

%     number of iterations for optimization fit
    iterFit = 3;
%     number of iterations for fit with random restart
    iterRandRestart = 0;


    fopts = optimoptions(fitFunc, 'Display', dispOpt, 'FunctionTolerance',...
                         funcTol, 'StepTolerance', stepTol,...
                         'OptimalityTolerance', optTol);


    varargout{1} = fitTransParams;
end

end
