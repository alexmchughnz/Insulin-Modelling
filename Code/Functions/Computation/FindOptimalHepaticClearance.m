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

DP = DebugPlots().FindOptimalHepaticClearance;

GRIDDEFAULTS = {[-0.1 0.775], [0.075 0.95], 0.025};

%% Setup
% Load grid settings.
if isempty(varargin)
    settings = GRIDDEFAULTS;
else
    settings = varargin;
end
nLBounds = settings{1};
xLBounds = settings{2};
delta = settings{3};

% Set up grid.
nLDelta = delta(1);
nLRange = nLBounds(1) : nLDelta : nLBounds(end);

xLDelta = delta(end);
xLRange = xLBounds(1) : xLDelta : xLBounds(end);

[xLGrid, nLGrid] = meshgrid(xLRange, nLRange);

if makeNewGrid || ~HasPersistent(P, "OptimalHepaticGrids")
    % Generate grid if we don't have one saved.  
    P = EvaluateGrid(P, nLGrid, xLGrid);
end
gridData = P.persistents.OptimalHepaticGrids{end};    


%% Find Optimal nL/xL
% Get minimum.
IResiduals = gridData.IResiduals;
nLGrid = gridData.nLGrid;
xLGrid = gridData.xLGrid;

[minIResidual, iiOptimal] = min(IResiduals(:));

P.results.nL = nLGrid(iiOptimal);
P.results.xL = xLGrid(iiOptimal);

% Get size of minimal error region.
deltaMSE = abs(IResiduals - minIResidual);
isWithin1SD = (deltaMSE <= P.persistents.stddevMSE);

optimalnL = nLGrid(isWithin1SD);
optimalxL = xLGrid(isWithin1SD);

[L, H] = bounds(optimalnL);
P.results.OptimalHepaticClearance.nLRange = [L H];

[L, H] = bounds(optimalxL);
P.results.OptimalHepaticClearance.xLRange = [L H];

P.results.OptimalHepaticClearance.minimalErrorRegionSize = sum(isWithin1SD(:));
P.results.OptimalHepaticClearance.minGridMSE = minIResidual;
%% ------------------------------------------------------------------------

%% Debug Plots
% Error Surface
if DP.ErrorSurface
    figTitle = sprintf("Error Surface @ %.3f", delta);
    MakeDebugPlot(figTitle, P, DP);
    
    % > Surface
    % Define surface parameters.
    nLRange = sort(unique(nLGrid));
    xLRange = sort(unique(xLGrid));
    
    gridMin = min(IResiduals(:));
    
    if isfield(P.data, 'stddevMSE')
        stddevMSE = P.data.stddevMSE;
    else
        stddevMSE = 20/100 * gridMin;  % Default to 20%.
    end
    
    % Define colors for different error regions.
    deltaMSE = abs(IResiduals - gridMin);
    isWithin1SD = (deltaMSE <= stddevMSE);
    isWithin3SD = (deltaMSE <= 3*stddevMSE);
    
    CO(:,:,1) = ~isWithin1SD * 1; %  red
    CO(:,:,2) = ones(size(deltaMSE)) * 1; % green
    CO(:,:,3) = ~isWithin3SD * 1; % blue
    
    % Plot surface.
    surf(xLRange, nLRange, IResiduals, CO,...
        'HandleVisibility', 'off', ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp');
    
    % > Contour
    % Define contour.
    numLevels = 10;
    minError = min(IResiduals(:));
%     levels = logspace(log10(min(IResiduals(:))), log10(max(IResiduals(:))), numLevels); % non-linear spacing
    levels = logspace(log10(minError), log10(1e+3*minError), numLevels); % non-linear spacing
    %             levels = linspace(min(IResiduals(:)), max(IResiduals(:)), numLevels); % linear spacing
    
    % Plot contour.
    contour3(xLRange, nLRange, IResiduals, ...
        levels, ...
        'Color', 'r', ...
        'HandleVisibility', 'off');
    
    
    % > Prettying
    delta = nLRange(2) - nLRange(1);
    
    xlim([0.1 0.9])
    ylim([0 0.7])
    
    xlabel("$x_L$ [min$^{-1}$]")
    ylabel("$n_L$")
    zlabel("Mean of squared errors [(mU/min)^2]") 
    
    % Redraw grid lines.
    spacing = 0.05;
    for ii = 1 : round(spacing/delta) : length(nLRange)
        plt = plot3(xLRange, ones(size(xLRange)) * nLRange(ii), IResiduals(ii,:));
        plt.Color = [0 0 0 0.1];
        plt.LineWidth = 0.2;
    end
    for ii = 1 : round(spacing/delta) : length(xLRange)
        plt = plot3(ones(size(nLRange)) * xLRange(ii), nLRange, IResiduals(:,ii));
        plt.Color = [0 0 0 0.1];
        plt.LineWidth = 0.2;
    end
    
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

%% Functions
function P = EvaluateGrid(P, nLGrid, xLGrid)

ISimulated = zeros([size(nLGrid) length(P.results.tArray)]);
IResiduals = zeros(size(nLGrid));

runtime = tic;
for ii = 1:numel(nLGrid)    
    message = sprintf('Searching at nL/xL = %g/%g...', nLGrid(ii), xLGrid(ii));
    PrintStatusUpdate(P, message, true);
    
    % Apply nL/xL for iteration.
    P.results.nL = nLGrid(ii);
    P.results.xL = xLGrid(ii);
    
    % Get other parameters and forward simulate models.
    P = FindGutEmptyingRate(P);
    P = FitInsulinSensitivity(P, false);
    P = SolveSystem(P);
    
    % Determine error at raw data points only.
    if (P.source == "Detemir")
        [tI, vI] = GetSimTime(P, P.data.ITotal);  % Data [mU/L]
        simI = P.results.I + P.results.IDF;       % Sim [mU/L]
    else
        [tI, vI] = GetSimTime(P, P.data.I);  % Data [mU/L]
        simI = P.results.I;                      % Sim [mU/L]
    end      
  
    inSimTime = GetTimeIndex(tI, P.results.tArray);
    
    dt = diff(tI);
    
    cumIntISim = cumtrapz(tI, simI(inSimTime));
    intISim = (cumIntISim(2:end) - cumIntISim(1:end-1)) ./ dt;
    
    cumIntIData = cumtrapz(tI, vI);
    intIData = (cumIntIData(2:end) - cumIntIData(1:end-1)) ./ dt;    
    
    intIErrors = intISim - intIData;
    
    % Save residuals.
    IResiduals(ii) = sum(intIErrors.^2)/numel(vI);  % Mean Squared Errors
    [row, col] = ind2sub(size(ISimulated), ii);
    ISimulated(row, col, :) = simI(:);
    
    runtime = PrintTimeRemaining("FindOptimalHepaticClearance", ...
        runtime, ii, numel(nLGrid), P);
end

% Export results.
saveStruct = struct(...
    'nLGrid', nLGrid, ...
    'xLGrid', xLGrid, ...
    'IResiduals', IResiduals, ...
    'ISimulated', ISimulated);

if ~HasPersistent(P, "OptimalHepaticGrids")
    P.persistents.OptimalHepaticGrids = {};
end

P.persistents.OptimalHepaticGrids{end+1} = saveStruct;   

end
