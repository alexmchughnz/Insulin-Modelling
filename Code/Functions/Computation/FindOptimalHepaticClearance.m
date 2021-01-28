function P = FindOptimalHepaticClearance(P, method, varargin)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'grid' mode - very slow!
% INPUT:
%   P        - patient struct
%   method   - 'grid' to perform grid search
%              'load' to load previously-generated residuals data
%   varargin - with 'grid', {1} nL boundary, and
%                           {2} xL boundary to search over
%                           {3} desired [nL, xL] grid precision
%            - otherwise, the filename to load (optional)
% OUTPUT:
%   P   - modified patient struct with nL and xL

global DEBUGPLOTS
global FILEFORMAT

GRIDFORMAT = "grid nL[%g %g]@%g xL[%g %g]@%g";
FILEFORMAT = '%s_%s.mat';

GRIDDEFAULTS = {[-0.1 0.775], [0.075 0.95], 0.025};

%% Setup
if method == "grid"
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
    
    
    filename = sprintf(GRIDFORMAT, ...
        nLBounds, nLDelta, xLBounds, xLDelta);
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, filename);
    
    
else
    if ~isempty(varargin)
        % Load by name.
        filename = varargin{1};
        load(ResultsPath(sprintf(FILEFORMAT, P.patientCode, filename)), ...
            'nLGrid', 'xLGrid', 'IResiduals');
    else
        % Load highest resolution grid file we can find.
        files = dir(ResultsPath(sprintf("%s_grid *", P.patientCode)));        
        resolutions = [];
        for ii = 1:length(files)
            [~, name, ~] = fileparts(files(ii).name);
            stringbits = split(name, "@");
            resolutions(ii) = str2double(stringbits(end));
        end
        
        [~, iiBest] = min(resolutions);
        file = files(iiBest);
        filename = file.name;
        load(fullfile(file.folder, filename), ...
            'nLGrid', 'xLGrid', 'IResiduals');
    end
end

if method == "refine"
    boundary = 0.1;
    delta = 0.005;
    
    % Find range surrounding 1SD area.
    gridMin = min(IResiduals(:));
    deltaMSE = abs(IResiduals - gridMin);
    isWithin1SD = (deltaMSE <= P.data.stddevMSE);
    
    nLMin = RoundToMultiple(min(nLGrid(isWithin1SD)), boundary) - boundary;
    nLMax = RoundToMultiple(max(nLGrid(isWithin1SD)), boundary) + boundary;
    xLMin = RoundToMultiple(min(xLGrid(isWithin1SD)), boundary) - boundary;
    xLMax = RoundToMultiple(max(xLGrid(isWithin1SD)), boundary) + boundary;
    
    % Set up grid.
    nLRange = nLMin : delta : nLMax;
    xLRange = xLMin : delta : xLMax;
    
    [xLGrid, nLGrid] = meshgrid(xLRange, nLRange);
    
    filename = sprintf(GRIDFORMAT, ...
        [nLMin nLMax], delta, [xLMin xLMax], delta);
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, filename);
end

%% Find Optimal nL/xL
% Get minimum.
[minIResidual, iiOptimal] = min(IResiduals(:));
bestnL = nLGrid(iiOptimal);
bestxL = xLGrid(iiOptimal);

P.results.nL = bestnL * ones(size(P.results.tArray));
P.results.xL = bestxL * ones(size(P.results.tArray));
P.results.minGridMSE = minIResidual;

% Get size of minimal error region.
gridMin = min(IResiduals(:));
deltaMSE = abs(IResiduals - gridMin);
isWithin1SD = (deltaMSE <= P.data.stddevMSE);


optimalnL = nLGrid(isWithin1SD);
optimalxL = xLGrid(isWithin1SD);

[L, H] = bounds(optimalnL);
P.results.optimalnLRange = [L H];
[L, H] = bounds(optimalxL);
P.results.optimalxLRange = [L H];

P.results.minimalErrorRegionSize = sum(isWithin1SD(:));

%% ------------------------------------------------------------------------

%% Debug Plots
DP = DEBUGPLOTS.FindOptimalHepaticClearance;

% Error Surface
if DP.ErrorSurface
    figTitle = sprintf("Error Surface @ %.3f", delta);
    MakeDebugPlot(figTitle, P, DP);
    
    % > Surface
    % Define surface parameters.
    nLRange = sort(unique(nLGrid));
    xLRange = sort(unique(xLGrid));
    
    gridMin = min(IResiduals(:));
    gridMax = max(IResiduals(:));
    
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
function [IResiduals, simI] = EvaluateGrid(PArray, nLGrid, xLGrid, savename)

global FILEFORMAT

ISimulated = cell(size(nLGrid));
IResiduals = zeros(size(nLGrid));
for ii = 1:numel(nLGrid)
    if length(PArray) == 1
        copyP = PArray(1);
    else
        copyP = PArray(ii);
    end
    
    message = sprintf('Searching at nL/xL = %g/%g...\n', nLGrid(ii), xLGrid(ii));
    PrintStatusUpdate(P, message);
    
    % Apply nL/xL for iteration.
    copyP.results.nL = nLGrid(ii) * ones(size(copyP.results.tArray));
    copyP.results.xL = xLGrid(ii) * ones(size(copyP.results.tArray));
    
    % Get other parameters and forward simulate models.
    copyP = FindGutEmptyingRate(copyP);
    copyP = FitInsulinSensitivity(copyP, false);
    copyP = SolveSystem(copyP);
    
    % Determine error at raw data points only.
    if (copyP.source == "Detemir")
        [tI, vI] = GetSimTime(copyP, copyP.data.ITotal);  % Data [mU/L]
        simI = copyP.results.I + copyP.results.IDF;       % Sim [mU/L]
    else
        [tI, vI] = GetSimTime(copyP, copyP.data.I);  % Data [mU/L]
        simI = copyP.results.I;                      % Sim [mU/L]
    end  
    inSimTime = GetTimeIndex(tI, copyP.results.tArray);
    IErrors = simI(inSimTime) - vI;
    
    % Save residuals.
    IResiduals(ii) = sum(IErrors.^2)/numel(vI);  % Mean Squared Errors
    ISimulated{ii} = simI;
    
    EstimateTimeRemaining(ii, numel(nLGrid))
end

% Export results.
save(ResultsPath(sprintf(FILEFORMAT, PArray.patientCode, savename)), ...
    'nLGrid', 'xLGrid', 'IResiduals', 'ISimulated')

end
