function P = GridSearchParameters(P, integralSystemFunc, makeNewGrid, gridOptions)
% Finds optimal parameters using grid search.
% Runs a LOT of forward simulations - very slow!
% INPUT:
%   P           - patient struct
%   makeNewGrid - bool for if to load or re-simulate
%   gridOptions - struct with ranges and stepsizes
% OUTPUT:
%   P   - modified patient struct with optimal parameters.

GRIDDEFAULTS.range = {[0 1], [0 1]};
GRIDDEFAULTS.step = [0.1 0.025];


%% Setup
% Run integral system to get names.
[~, ~, ~, ~, paramNames] = integralSystemFunc(P);
gridName = paramNames(1) + paramNames(2) + "Grids";

if ~exist("gridOptions", "var")
    gridOptions = GRIDDEFAULTS;
end

[P, hasGrids] = GetPersistent(P, gridName);
if makeNewGrid || ~hasGrids
    % Load grid settings.
    param1Range = gridOptions.range{1};
    param1Step = gridOptions.step(1);
    param2Range = gridOptions.range{end};
    param2Step = gridOptions.step(end);
    
    % Set up grid.
    param1Range = param1Range(1) : param1Step : param1Range(end);
    param2Range = param2Range(1) : param2Step : param2Range(end);
    
    [param2Grid, param1Grid] = meshgrid(param2Range, param1Range);
    
    % Generate grid if we don't have one saved.
    P = EvaluateGrid(P, param1Grid, param2Grid, integralSystemFunc, gridName);
end

gridData = P.persistents.(gridName){end};


%% Find Optimal nL/param2
% Get minimum.
objectiveValues = gridData.objectiveValues;
param1Grid = gridData.param1Grid;
param2Grid = gridData.param2Grid;

[objectiveMin, iiMinimum] = min(objectiveValues(:));

P.results.(paramNames(1)) = param1Grid(iiMinimum);
P.results.(paramNames(2)) = param2Grid(iiMinimum);

% Get size of minimal error region.
deltaMSE = abs(objectiveValues - objectiveMin);
isWithin1SD = (deltaMSE <= P.persistents.stddevMSE);

optimalParam1 = param1Grid(isWithin1SD);
optimalParam2 = param2Grid(isWithin1SD);

[L, H] = bounds(optimalParam1);
P.results.GridSearch.(paramNames(1)+"Range") = [L H];

[L, H] = bounds(optimalParam2);
P.results.GridSearch.(paramNames(2)+"Range") = [L H];

P.results.GridSearch.minimalErrorRegionSize = sum(isWithin1SD(:));
P.results.GridSearch.minGridMSE = objectiveMin;


%% Plotting
plotvars.gridName = gridName;
plotvars.paramNames = paramNames;

MakePlots(P, plotvars);

end


function P = EvaluateGrid(P, param1Grid, param2Grid, integralSystemFunc, gridName)
runtime = tic;

[A, b, ~, ~, paramNames] = integralSystemFunc(P);

% Set up results grids.
objectiveValues = zeros(size(param1Grid));
integralErrors = zeros(size(param1Grid));
dataErrors = zeros(size(param1Grid));

% Get integral system for patient.
for ii = 1:numel(param1Grid)
    message = sprintf('Searching at %s/%s = %g/%g...', paramNames(1), paramNames(2), param1Grid(ii), param2Grid(ii));
    PrintStatusUpdate(P, message, true);
    
    % Apply param1/param2 for iteration.
    P.results.(paramNames(1)) = param1Grid(ii);
    P.results.(paramNames(2)) = param2Grid(ii);
    
    % Get other parameters and forward simulate models.
    P = FindGutEmptyingRate(P);
    P = FitInsulinSensitivity(P);
    P = SolveSystem(P);
    
    % Determine error.
    x = [P.results.(paramNames(1));
         P.results.(paramNames(2))];
    integralError = sum((A*x - b).^2);
    
    [tI, vI] = GetData(P.data.I);
    [~, simI] = GetResultsSample(P, tI, P.results.I);
    dataError = sum((vI-simI).^2);
    
    totalObjectiveValue = integralError;
    
    % Save residuals.
    integralErrors(ii) = integralError;
    dataErrors(ii) = dataError;
    objectiveValues(ii) = totalObjectiveValue;
    
    runtime = PrintTimeRemaining("GridSearchParameters", ...
        runtime, ii, numel(param1Grid), P);
end

% Export results.
saveStruct = struct(...
    'param1Grid', param1Grid, ...
    'param2Grid', param2Grid, ...
    'integralErrors', integralErrors, ...
    'dataErrors', dataErrors, ...
    'objectiveValues', objectiveValues);

[P, hasGrids] = GetPersistent(P, gridName);
if ~hasGrids
    P.persistents.(gridName) = {};
end

P.persistents.(gridName){end+1} = saveStruct;

end


function MakePlots(P, plotvars)
DP = DebugPlots().GridSearchParameters;

gridData = P.persistents.(plotvars.gridName){end};
objectiveValues = gridData.objectiveValues;
param1Grid = gridData.param1Grid;
param2Grid = gridData.param2Grid;


%% Error Surface
if DP.ErrorSurface
    figTitle = sprintf("%s-%s Error Surface", plotvars.paramNames(1), plotvars.paramNames(2));
    MakeDebugPlot(figTitle, P, DP);
    
    %% Surface
    % Define surface parameters.
    param1Range = sort(unique(param1Grid));
    param2Range = sort(unique(param2Grid));
    gridMin = min(objectiveValues(:));
    
    if isfield(P.data, 'stddevMSE')
        stddevMSE = P.data.stddevMSE;
    else
        stddevMSE = 20/100 * gridMin;  % Default to 20%.
    end
    
    % Define colors for different error regions.
    deltaMSE = abs(objectiveValues - gridMin);
    isWithin1SD = (deltaMSE <= stddevMSE);
    isWithin3SD = (deltaMSE <= 3*stddevMSE);
    
    CO(:,:,1) = ~isWithin1SD * 1; %  red
    CO(:,:,2) = ones(size(deltaMSE)) * 1; % green
    CO(:,:,3) = ~isWithin3SD * 1; % blue
    
    % Plot surface.
    surf(param2Range, param1Range, objectiveValues, CO,...
        'HandleVisibility', 'off', ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');
    
    %% Contour
    % Define contour.
    numLevels = 10;
    minError = min(objectiveValues(:));
    levels = logspace(log10(minError), log10(1e+3*minError), numLevels); % non-linear spacing
    
    % Plot contour.
    contour3(param2Range, param1Range, objectiveValues, ...
        levels, ...
        'Color', 'r', ...
        'HandleVisibility', 'off');
    
    %% Prettying
    xlim([0 1])
    ylim([0 1])
    
    xlabel(plotvars.paramNames(2))
    ylabel(plotvars.paramNames(1))
    zlabel("Mean of squared errors [(mU/min)^2]")
    
    % Add physiological region.
    param2Phys = [0.5 0.9];
    param1Phys = [0.1 0.3];
    x = param2Phys([1 1 end end]);
    y = param1Phys([1 end end 1]);
    z = 1e+6 * ones(1, 4);
    patch(x, y, z, 'r',...
        'FaceColor', '#D95319', ...
        'FaceAlpha', 0.2, ...
        'EdgeColor', 'none')
end
end

