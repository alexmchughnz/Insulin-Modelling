function P = FindOptimalHepaticClearance(P, method, varargin)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'grid' mode - very slow!
% INPUT:
%   P        - patient struct
%   method   - 'grid' to perform grid search
%              'load' to load previously-generated residuals data
%   varargin - with 'grid', {1} nL bounds, and
%                           {2} xL bounds to search over
%                           {3} desired [nL, xL] grid precision
%            - with 'load', the filename to load (optional)
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
    bounds = 0.05;
    delta = 0.005;
    
    % Find range surrounding 1SD area.
    gridMin = min(IResiduals(:));
    deltaMSE = abs(IResiduals - gridMin);
    isWithin1SD = (deltaMSE <= P.data.stddevMSE);
    
    nLMin = RoundToMultiple(min(nLGrid(isWithin1SD)), bounds) - bounds;
    nLMax = RoundToMultiple(max(nLGrid(isWithin1SD)), bounds) + bounds;
    xLMin = RoundToMultiple(min(xLGrid(isWithin1SD)), bounds) - bounds;
    xLMax = RoundToMultiple(max(xLGrid(isWithin1SD)), bounds) + bounds;
    
    % Set up grid.
    nLRange = nLMin : delta : nLMax;
    xLRange = xLMin : delta : xLMax;
    
    [xLGrid, nLGrid] = meshgrid(xLRange, nLRange);
    
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, filename);
end

%% Find Optimal nL/xL
[minIResidual, iiOptimal] = min(IResiduals(:));
bestnL = nLGrid(iiOptimal);
bestxL = xLGrid(iiOptimal);

P.results.nL = bestnL * ones(size(P.results.tArray));
P.results.xL = bestxL * ones(size(P.results.tArray));
P.results.minGridMSE = minIResidual;

%% ------------------------------------------------------------------------

%% Debug Plots
DP = DEBUGPLOTS.FindOptimalHepaticClearance;

% Error Surface
if DP.ErrorSurface
    MakeDebugPlot(P, DP);
    hold on
    
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
    
    CO(:,:,1) = (isWithin3SD & ~isWithin1SD) * 0.5; %  red
    CO(:,:,2) = (isWithin3SD) * 0.5; % green
    CO(:,:,3) = ~(isWithin3SD) .* (gridMax-IResiduals)/gridMax; % blue
    
    % Plot surface.
    surf(xLRange, nLRange, IResiduals, CO,...
        'HandleVisibility', 'off', ...
        'FaceColor', 'interp');
    
    % > Contour
    % Define contour.
    numLevels = 15;
    levels = logspace(log10(min(IResiduals(:))), log10(max(IResiduals(:))), numLevels); % non-linear spacing
    %             levels = linspace(min(IResiduals(:)), max(IResiduals(:)), numLevels); % linear spacing
    
    % Plot contour.
    contour3(xLRange, nLRange, IResiduals, ...
        levels, ...
        'Color', 'r', ...
        'HandleVisibility', 'off');
    
    title(sprintf("%s: Error Surface", P.patientCode))
    
    xlabel("$x_L$ [1/min]")
    ylabel("$n_L$ [-]")
    zlabel("Mean of squared errors [(mU/min)^2]")
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
    
    fprintf('\nP%d: Searching at nL/xL = %g/%g...\n', ...
        copyP.patientNum, nLGrid(ii), xLGrid(ii))
    
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
