function P = FindOptimalHepaticClearance(P, makeNewGrid, varargin)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'grid' mode - very slow!
% INPUT:
%   P        - patient struct
%   varargin - {1} nL boundary, and
%              {2} xL boundary to search over
%              {3} desired [nL, xL] grid precision
% OUTPUT:
%   P   - modified patient struct with nL and xL

GRIDDEFAULTS = {[-0.1 0.775], [0.075 0.95], 0.02};

%% Setup
if makeNewGrid || ~HasPersistent(P, "OptimalHepaticGrids")
    % Load grid settings.
    if isempty(varargin)
        settings = GRIDDEFAULTS;
    else
        settings = varargin;
    end
    nLBounds = settings{1};
    xLBounds = settings{2};
    nLxLDelta = settings{3};
    
    % Set up grid.
    nLDelta = nLxLDelta(1);
    nLRange = nLBounds(1) : nLDelta : nLBounds(end);
    
    xLDelta = nLxLDelta(end);
    xLRange = xLBounds(1) : xLDelta : xLBounds(end);
    
    [xLGrid, nLGrid] = meshgrid(xLRange, nLRange);
    
    % Generate grid if we don't have one saved.
    P = EvaluateGrid(P, nLGrid, xLGrid);
end
gridData = P.persistents.OptimalHepaticGrids{end};


%% Find Optimal nL/xL
% Get minimum.
objectiveValues = gridData.objectiveValues;
nLGrid = gridData.nLGrid;
xLGrid = gridData.xLGrid;

[objectiveMin, iiMinimum] = min(objectiveValues(:));

P.results.nL = nLGrid(iiMinimum);
P.results.xL = xLGrid(iiMinimum);

% Get size of minimal error region.
deltaMSE = abs(objectiveValues - objectiveMin);
isWithin1SD = (deltaMSE <= P.persistents.stddevMSE);

optimalnL = nLGrid(isWithin1SD);
optimalxL = xLGrid(isWithin1SD);

[L, H] = bounds(optimalnL);
P.results.OptimalHepaticClearance.nLRange = [L H];

[L, H] = bounds(optimalxL);
P.results.OptimalHepaticClearance.xLRange = [L H];

P.results.OptimalHepaticClearance.minimalErrorRegionSize = sum(isWithin1SD(:));
P.results.OptimalHepaticClearance.minGridMSE = objectiveMin;


%% Plotting
MakePlots(P);

% Error Comparison
if DP.ErrorComparison
    figTitle = "Error Comparison";
    MakeDebugPlot(figTitle, P, DP);
    
    % Integral Error
    subplot(1,2,1)
    hold on
    
    % Define surface parameters.
    nLRange = sort(unique(nLGrid));
    xLRange = sort(unique(xLGrid));  
    
    % Plot surface.
    surf(xLRange, nLRange, gridData.integralErrors, ...
        'HandleVisibility', 'off', ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');
    
    % Plot contour.
    minError = min(gridData.integralErrors(:));
    levels = logspace(log10(minError), log10(1e+3*minError), numLevels); % non-linear spacing
    contour3(xLRange, nLRange, gridData.integralErrors, ...
        levels, ...
        'Color', 'r', ...
        'HandleVisibility', 'off');
    
    xlim([0 1])
    ylim([0 1])
    
    
    % Data Error
    subplot(1,2,2)
    hold on
    
    % Define surface parameters.
    nLRange = sort(unique(nLGrid));
    xLRange = sort(unique(xLGrid));  
    
    % Plot surface.
    surf(xLRange, nLRange, gridData.dataErrors, ...
        'HandleVisibility', 'off', ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');    
    
    % Plot contour.
    minError = min(gridData.dataErrors(:));
    levels = logspace(log10(minError), log10(1e+3*minError), numLevels);
    contour3(xLRange, nLRange, gridData.dataErrors, ...
        levels, ...
        'Color', 'r', ...
        'HandleVisibility', 'off');
    
    xlim([0 1])
    ylim([0 1])
end


end


function P = EvaluateGrid(P, nLGrid, xLGrid)
runtime = tic;

ISimulated = zeros([size(nLGrid) length(P.results.tArray)]);
objectiveValues = zeros(size(nLGrid));
integralErrors = zeros(size(nLGrid));
dataErrors = zeros(size(nLGrid));

% Get integral system for patient.
A = P.results.integrals.A;
b = P.results.integrals.b;

for ii = 1:numel(nLGrid)
    message = sprintf('Searching at nL/xL = %g/%g...', nLGrid(ii), xLGrid(ii));
    PrintStatusUpdate(P, message, true);
    
    % Apply nL/xL for iteration.
    P.results.nL = nLGrid(ii);
    P.results.xL = xLGrid(ii);
    
    % Get other parameters and forward simulate models.
    P = FindGutEmptyingRate(P);
    P = FitInsulinSensitivity(P);
    P = SolveSystem(P);
    
    % Determine error.
    x = [P.results.nL; 1 - P.results.xL];
    integralError = sum((A*x - b).^2);
    
    [tI, vI] = GetData(P.data.I);
    [~, simI] = GetResultsSample(P, tI, P.results.I);
    dataError = sum((vI-simI).^2);
    
    scale = 0;
    totalObjectiveValue = integralError + scale*dataError;
    
    % Save residuals.
    integralErrors(ii) = integralError;
    dataErrors(ii) = dataError;
    objectiveValues(ii) = totalObjectiveValue;
    [row, col] = ind2sub(size(P.results.I), ii);
    ISimulated(row, col, :) = P.results.I(:);
    
    runtime = PrintTimeRemaining("FindOptimalHepaticClearance", ...
        runtime, ii, numel(nLGrid), P);
end

% Export results.
saveStruct = struct(...
    'nLGrid', nLGrid, ...
    'xLGrid', xLGrid, ...
    'integralErrors', integralErrors, ...
    'dataErrors', dataErrors, ...
    'objectiveValues', objectiveValues, ...
    'ISimulated', ISimulated);

if ~HasPersistent(P, "OptimalHepaticGrids")
    P.persistents.OptimalHepaticGrids = {};
end

P.persistents.OptimalHepaticGrids{end+1} = saveStruct;

end


function MakePlots(P)
DP = DebugPlots().FindOptimalHepaticClearance;

gridData = P.persistents.OptimalHepaticGrids{end};
objectiveValues = gridData.objectiveValues;
nLGrid = gridData.nLGrid;
xLGrid = gridData.xLGrid;


%% Error Surface
if DP.ErrorSurface
    figTitle = "Error Surface";
    MakeDebugPlot(figTitle, P, DP);
    
    %% Surface
    % Define surface parameters.
    nLRange = sort(unique(nLGrid));
    xLRange = sort(unique(xLGrid));
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
    surf(xLRange, nLRange, objectiveValues, CO,...
        'HandleVisibility', 'off', ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');
    
    %% Contour
    % Define contour.
    numLevels = 10;
    minError = min(objectiveValues(:));
    levels = logspace(log10(minError), log10(1e+3*minError), numLevels); % non-linear spacing
    
    % Plot contour.
    contour3(xLRange, nLRange, objectiveValues, ...
        levels, ...
        'Color', 'r', ...
        'HandleVisibility', 'off');
    
    %% Prettying
    xlim([0 1])
    ylim([0 1])
    
    xlabel("$x_L$ [min$^{-1}$]")
    ylabel("$n_L$")
    zlabel("Mean of squared errors [(mU/min)^2]")
    
    % Add physiological region.
    xLPhys = [0.5 0.9];
    nLPhys = [0.1 0.3];
    x = xLPhys([1 1 end end]);
    y = nLPhys([1 end end 1]);
    z = 1e+6 * ones(1, 4);
    patch(x, y, z, 'r',...
        'FaceColor', '#D95319', ...
        'FaceAlpha', 0.2, ...
        'EdgeColor', 'none')
end
end

