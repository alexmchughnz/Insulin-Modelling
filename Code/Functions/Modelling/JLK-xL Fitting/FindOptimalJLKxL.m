function P = FindOptimalJLKxL(P, makeNewGrid)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'grid' mode - very slow!
% INPUT:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with nL and xL

GRIDNAME = "JLKxLGrids";
GRIDDEFAULTS = {[0.1 1.0], [0.075 0.95], [0.05 0.025]};

%% Setup
[P, hasGrids] = GetPersistent(P, GRIDNAME);
if makeNewGrid || ~hasGrids
    % Load grid settings.
    settings = GRIDDEFAULTS;
    JLKBounds = settings{1};
    xLBounds = settings{2};
    JLKxLDelta = settings{3};
    
    % Set up grid.
    JLKDelta = JLKxLDelta(1);
    JLKRange = JLKBounds(1) : JLKDelta : JLKBounds(end);
    
    xLDelta = JLKxLDelta(end);
    xLRange = xLBounds(1) : xLDelta : xLBounds(end);
    
    [xLGrid, JLKGrid] = meshgrid(xLRange, JLKRange);
    
    % Generate grid if we don't have one saved.
    P = EvaluateGrid(P, JLKGrid, xLGrid, GRIDNAME);
end
gridData = P.persistents.(GRIDNAME){end};


%% Find Optimal nL/xL
% Get minimum.
objectiveValues = gridData.objectiveValues;
JLKGrid = gridData.JLKGrid;
xLGrid = gridData.xLGrid;

[objectiveMin, iiMinimum] = min(objectiveValues(:));

P.results.JLK = JLKGrid(iiMinimum);
P.results.xL = xLGrid(iiMinimum);

% Get size of minimal error region.
deltaMSE = abs(objectiveValues - objectiveMin);
isWithin1SD = (deltaMSE <= P.persistents.stddevMSE);

optimalJLK = JLKGrid(isWithin1SD);
optimalxL = xLGrid(isWithin1SD);

[L, H] = bounds(optimalJLK);
P.results.FindOptimal.JLKRange = [L H];

[L, H] = bounds(optimalxL);
P.results.FindOptimal.xLRange = [L H];

P.results.FindOptimal.minimalErrorRegionSize = sum(isWithin1SD(:));
P.results.FindOptimal.minGridMSE = objectiveMin;


%% Plotting
plotvars.GRIDNAME = GRIDNAME;

MakePlots(P, plotvars);

end


function P = EvaluateGrid(P, JLKGrid, xLGrid, GRIDNAME)
runtime = tic;

ISimulated = zeros([size(JLKGrid) length(P.results.tArray)]);
objectiveValues = zeros(size(JLKGrid));
integralErrors = zeros(size(JLKGrid));
dataErrors = zeros(size(JLKGrid));

% Get integral system for patient.
A = P.results.integrals.A;
b = P.results.integrals.b;

for ii = 1:numel(JLKGrid)
    message = sprintf('Searching at JLK/xL = %g/%g...', JLKGrid(ii), xLGrid(ii));
    PrintStatusUpdate(P, message, true);
    
    % Apply JLK/xL for iteration.
    P.results.JLK = JLKGrid(ii);
    P.results.xL = xLGrid(ii);
    
    % Get other parameters and forward simulate models.
    P = FindGutEmptyingRate(P);
    P = FitInsulinSensitivity(P);
    P = SolveSystem(P);
    
    % Determine error.
    x = [P.results.JLK; 1 - P.results.xL];
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
    
    runtime = PrintTimeRemaining("FindFindOptimal", ...
        runtime, ii, numel(JLKGrid), P);
end

% Export results.
saveStruct = struct(...
    'JLKGrid', JLKGrid, ...
    'xLGrid', xLGrid, ...
    'integralErrors', integralErrors, ...
    'dataErrors', dataErrors, ...
    'objectiveValues', objectiveValues, ...
    'ISimulated', ISimulated);

[P, hasGrids] = GetPersistent(P, GRIDNAME);
if ~hasGrids
    P.persistents.(GRIDNAME) = {};
end

P.persistents.(GRIDNAME){end+1} = saveStruct;

end


function MakePlots(P, plotvars)
DP = DebugPlots().FindOptimal;

gridData = P.persistents.(plotvars.GRIDNAME){end};
objectiveValues = gridData.objectiveValues;
JLKGrid = gridData.JLKGrid;
xLGrid = gridData.xLGrid;


%% Error Surface
if DP.ErrorSurface
    figTitle = "JLK-xL Error Surface";
    MakeDebugPlot(figTitle, P, DP);
    
    %% Surface
    % Define surface parameters.
    JLKRange = sort(unique(JLKGrid));
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
    surf(xLRange, JLKRange, objectiveValues, CO,...
        'HandleVisibility', 'off', ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');
    
    %% Contour
    % Define contour.
    numLevels = 10;
    minError = min(objectiveValues(:));
    levels = logspace(log10(minError), log10(1e+3*minError), numLevels); % non-linear spacing
    
    % Plot contour.
    contour3(xLRange, JLKRange, objectiveValues, ...
        levels, ...
        'Color', 'r', ...
        'HandleVisibility', 'off');
    
    %% Prettying
    xlim([0 1])
    ylim([0 1])
    
    xlabel("$x_L$ [min$^{-1}$]")
    ylabel("$J_LK$")
    zlabel("Mean of squared errors [(mU/min)^2]")
    
    % Add physiological region.
    xLPhys = [0.5 0.9];
    JLKPhys = [0.1 0.3];
    x = xLPhys([1 1 end end]);
    y = JLKPhys([1 end end 1]);
    z = 1e+6 * ones(1, 4);
    patch(x, y, z, 'r',...
        'FaceColor', '#D95319', ...
        'FaceAlpha', 0.2, ...
        'EdgeColor', 'none')
end
end

