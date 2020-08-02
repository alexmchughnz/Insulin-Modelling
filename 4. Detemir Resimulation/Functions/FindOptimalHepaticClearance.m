function P = FindOptimalHepaticClearance(P, method, varargin)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'grid' mode - very slow!
% INPUT:
%   P        - patient struct
%   method   - 'grid' to perform grid search
%              'line' to perform grid search on line
%              '2dline' to perform 2D search on line
%              'load' to load previously-generated residuals data
%              'improve' to load data and iterate it to some precision
%   varargin - with 'grid', {1} nL bounds, and
%                           {2} xL bounds to search over
%                           {3} desired [nL, xL] grid precision
%            - with 'line'/'2dline', {1} nL-intercept, and
%                                    {2} xL-intercept to search between
%                                    {3} desired [nL, xL] grid precision
%            - with 'load', the filename to load
%            - with 'improve', {1} the filename to load and improve
%                              {2} desired [nL, xL] grid precision
% OUTPUT:
%   P   - modified patient struct with nL and xL

global DEBUGPLOTS
global FILEFORMAT

GRIDFORMAT = "grid nL[%g %g]@%g xL[%g %g]@%g";
LINEFORMAT = "line nL=%g@%g to xL=%g@%g, t=%g";
LINE2DFORMAT = "2dline nL=%g@%g to xL=%g";
FILEFORMAT = '%s_%s.mat';

nLtoxLLineFun = @(nLIntercept, xLIntercept, nLDelta) ...
    (@(nL) RoundToMultiple(xLIntercept - xLIntercept/nLIntercept * nL, ...
    nLDelta));

%% Setup
if isequal(method, 'grid')
    type = 'grid';
    
    nLBounds = varargin{1};
    xLBounds = varargin{2};
    delta = varargin{3};
    
    nLDelta = delta(1);
    nLRange = nLBounds(1) : nLDelta : nLBounds(end);
    
    xLDelta = delta(end);
    xLRange = xLBounds(1) : xLDelta : xLBounds(end);
    
    [xLGrid, nLGrid] = meshgrid(xLRange, nLRange);
    
    savename = sprintf(GRIDFORMAT, ...
        nLBounds, nLDelta, xLBounds, xLDelta);
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename);
    
elseif isequal(method, 'line')
    type = 'line';
    
    nLIntercept = varargin{1};
    xLIntercept = varargin{2};
    delta = varargin{3};
    
    % Establish nL range.
    nLDelta = delta(1);
    xLDelta = delta(end);
    
    nLBounds = [0 1];
    nLRange = nLBounds(1) : nLDelta : nLBounds(end);
    
    % Define line - a series of slices of xL for each nL value.
    nLtoxL = nLtoxLLineFun(nLIntercept, xLIntercept, nLDelta);
    
    thickness = 0.2;  % In xL axis [-]
    xLLine = [];
    nLLine = [];
    for nL = nLRange
        xL = nLtoxL(nL);  % Corresponding xL on line for nL.
        xLSlice = xL-thickness/2 : xLDelta : xL+thickness/2;
        nLSlice = nL * ones(size(xLSlice));
        
        % Append to values along line.
        xLLine = [xLLine xLSlice];
        nLLine = [nLLine nLSlice];
    end
    
    % Find residuals along nL/xL ranges.
    savename = sprintf(LINEFORMAT, ...
        nLIntercept, nLDelta, xLIntercept, xLDelta, thickness);
    LineResiduals = EvaluateGrid(P, nLLine, xLLine, savename);
    
    % Reshape residuals onto grid.
    xLRange = sort(unique(xLLine));
    
    [xLGrid, nLGrid] = meshgrid(xLRange, nLRange);
    IResiduals = nan(size(nLGrid));
    for ii = 1 : length(LineResiduals)
        iinL = find(nLRange == nLLine(ii));
        iixL = find(xLRange == xLLine(ii));
        
        IResiduals(iixL, iinL) = LineResiduals(ii);
    end
    
    
elseif isequal(method, '2dline')
    type = '2dline';
    
    nLIntercept = varargin{1};
    xLIntercept = varargin{2};
    delta = varargin{3};
    
    % Establish full grid over area.
    nLDelta = delta(1);
    nLBounds = [0 1];
    nLGrid = nLBounds(1) : nLDelta : nLBounds(end);
    
    nLtoxL = nLtoxLLineFun(nLIntercept, xLIntercept, nLDelta);
    xLGrid = nLtoxL(nLGrid);
    
    % Find residuals along nL/xL line.
    savename = sprintf(LINE2DFORMAT, ...
        nLIntercept, nLDelta, xLIntercept);
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename);
    
elseif isequal(method, 'load')
    % Load by name.
    loadname = varargin{1};
    load(ResultsPath(sprintf(FILEFORMAT, P.patientCode, loadname)), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    
    words = split(loadname, ' ');
    type = words{1};
    if type == "line"
        % Reshape residuals onto grid.
        nLLine = nLGrid;
        nLRange = sort(unique(nLLine));
        xLLine = xLGrid;
        xLRange = sort(unique(xLLine));
        [xLGrid, nLGrid] = meshgrid(xLRange, nLRange);
        
        LineResiduals = IResiduals;
        
        IResiduals = nan(numel(nLRange), numel(xLRange));
        for ii = 1 : length(IResiduals)
            iinL = find(nLRange == nLLine(ii));
            iixL = find(xLRange == xLLine(ii));
            
            IResiduals(iinL, iixL) = LineResiduals(ii);
        end
        
    elseif type == "2dline"
        nLIntercept = nLGrid(xLGrid == 0);
        xLIntercept = xLGrid(nLGrid == 0);
        nLDelta = diff(nLGrid(1:2));
        
        nLtoxL = nLtoxLLineFun(nLIntercept, xLIntercept, nLDelta);
        
    else
        nLRange = sort(unique(nLGrid));
        xLRange = sort(unique(xLGrid));
    end
    
    
elseif isequal(method, 'improve')
    type = 'grid';
    
    % Load by name.
    loadname = varargin{1};
    delta = varargin{2};
    nLPrecision = delta(1);
    xLPrecision = delta(end);
    
    load(ResultsPath(sprintf(FILEFORMAT, P.patientCode, loadname)), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    
    nLRange = sort(unique(nLGrid));
    xLRange = sort(unique(xLGrid));
    nLDelta = diff(nLRange(1:2));
    xLDelta = diff(xLRange(1:2));
    
    loop = true;
    factor = 10;
    numSegments = 5;
    while loop
        % Find optima of current grid.
        iiOpt = find(IResiduals == min(IResiduals(:)));
        [iixLOpt, iinLOpt] = ind2sub(size(IResiduals), iiOpt);
        nLOpt = nLRange(iinLOpt);
        xLOpt = xLRange(iixLOpt);
        
        % Set up new search range around optima.
        nLBounds = [nLOpt - nLDelta,  nLOpt + nLDelta];
        xLBounds = [xLOpt - xLDelta,  xLOpt + xLDelta];
        nLRange = linspace(nLBounds(1), nLBounds(end), numSegments);
        xLRange = linspace(xLBounds(1), xLBounds(end), numSegments);
        [xLGrid, nLGrid] = meshgrid(xLRange, nLRange);
        
        % Evaluate over grid.
        savename = sprintf(GRIDFORMAT, ...
            nLBounds, nLDelta, xLBounds, xLDelta);
        IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename);
        
        % Update grid settings for next iteration.
        if (iinLOpt == 1) || (iinLOpt == numSegments) ...
                || (iixLOpt == 1) || (iixLOpt == numSegments)
            % Best value is on boundary. Next loop will recenter window
            % at same precision - do nothing and continue!
            loop = true;
            
        elseif (nLDelta == nLPrecision) && (xLDelta == xLPrecision)
            % We're at desired precision, and optimum isn't on edge.
            % Time to finish up.
            loop = false;
            
        else
            % Increase precision of grid.
            if nLDelta/factor > nLPrecision
                nLDelta = nLDelta/factor;
            else
                nLDelta = nLPrecision;
            end
            
            if xLDelta/factor > xLPrecision
                xLDelta = xLDelta/factor;
            else
                xLDelta = xLPrecision;
            end
        end
    end
    
    
elseif isequal(method, 'variance')
    % Load by name.
    loadname = varargin{1};
    variance = varargin{2};
    
    load(ResultsPath(sprintf(FILEFORMAT, P.patientCode, loadname)), ...
        'nLGrid', 'xLGrid', 'IResiduals', 'ISimulated');
    
    words = split(loadname, ' ');
    type = words{1};
    
    % Prepare array of patients for each nL/xL point.
    bestPArray = [];
    worstPArray = [];
    
    [~, vITotal] = GetSimTime(P, P.data.ITotal);  % Data [pmol/L]
    
    for ii = length(ISimulated)
        isModelOver = (ISimulated{ii} > vITotal);
        
        bestDataVariance = (1+variance)*isModelOver + (1-variance)*~isModelOver;
        worstDataVariance = (1-variance)*isModelOver + (1+variance)*~isModelOver;
        
        bestP = P;
        bestP.data.ITotal.value = bestP.data.ITotal.value .* bestDataVariance;
        bestPArray = [bestPArray bestP];
        
        worstP = P;
        worstP.data.ITotal.value = worstP.data.ITotal.value .* worstDataVariance;
        worstPArray = [worstPArray worstP];
    end
    
    % Evaluate original nL/xL range for both worst and best cases.
    savename = sprintf("%s best", loadname);  % Best Case
    EvaluateGrid(bestPArray, nLGrid, xLGrid, savename);
    
    savename = sprintf("%s worst", loadname);  % Worst Case
    EvaluateGrid(worstPArray, nLGrid, xLGrid, savename);
    
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
    if ~exist('loadname', 'var')
        loadname = savename;
    end
    
    if type == "2dline"
        % Load original.
        originalResiduals = IResiduals;
        plt = plot(nLGrid, originalResiduals);
        plt.DisplayName = 'Residuals';
        
        % Load best/worst files if they exist.
        bestFile = sprintf(FILEFORMAT, P.patientCode, loadname + " best");
        if exist(bestFile, 'file')
            load(ResultsPath(bestFile), ...
                'nLGrid', 'IResiduals');
            bestResiduals = IResiduals;
            plt = plot(nLGrid, bestResiduals);
            plt.DisplayName = 'Best Case Data Residuals';
        end
        
        worstFile = sprintf(FILEFORMAT, P.patientCode, loadname + " worst");
        if exist(worstFile, 'file')
            load(ResultsPath(worstFile), ...
                'nLGrid', 'IResiduals');
            worstResiduals = IResiduals;
            plt = plot(nLGrid, worstResiduals);
            plt.DisplayName = 'Worst Case Data Residuals';
        end
        
        plt = plot(nLGrid(iiOptimal), minIResidual, 'r*');
        plt.DisplayName = 'Optimal Point';
        
        legend
        
        ax1 = gca;
        ax1.XColor = 'b';
        ax1XTick = ax1.XTick';
        
        ax2 = axes();
        ax2.XAxisLocation = 'top';
        ax2.XDir = 'reverse';
        ax2.Color = 'none';
        ax2.XColor = 'r';
        ax2.XTick = linspace(0, 1, length(ax1XTick));
        ax2.XTickLabel = num2str(flip(nLtoxL(ax1XTick)));
        ax2.YAxis.Visible = 'off';
        
        figTitle = sprintf("%s: Error Line of (I+IDF) Fitting Along $x_L = %.2f - %.2f n_L$", ...
            P.patientCode, xLIntercept, xLIntercept/nLIntercept);
        title(figTitle)
        xlabel(ax1, "$n_L$ [-]");
        xlabel(ax2, "$x_L$ [1/min]");
        ylabel(ax1, "2-norm of residuals, $\psi$ [mU/min]")
        
        grid on
        ax1.Position = ax2.Position;
        
    else
        surfaces = {};
        
        % Queue original surface.
        S = struct('IResiduals', IResiduals, ...
            'nLGrid', nLGrid, ...
            'xLGrid', xLGrid, ...
            'name', 'Residuals');
        surfaces = [surfaces S];
        
        % Queue best/worst if they exist.
        bestFile = ResultsPath(sprintf(FILEFORMAT, P.patientCode, loadname + " best"));
        if exist(bestFile, 'file')
            load(ResultsPath(bestFile), ...
                'nLGrid', 'IResiduals');
            S = struct('IResiduals', IResiduals, ...
                'nLGrid', nLGrid, ...
                'xLGrid', xLGrid, ...
                'name', 'Best Case Data Residuals');
            surfaces = [surfaces S];
        end
        
        worstFile = ResultsPath(sprintf(FILEFORMAT, P.patientCode, loadname + " worst"));
        if exist(worstFile, 'file')
            load(ResultsPath(worstFile), ...
                'nLGrid', 'IResiduals');
            S = struct('IResiduals', IResiduals, ...
                'nLGrid', nLGrid, ...
                'xLGrid', xLGrid, ...
                'name', 'Worst Case Data Residuals');
            surfaces = [surfaces S];
        end
        
        % Plot each surface.
        for ii = 1:length(surfaces)
            nLRange = sort(unique(S.nLGrid));
            xLRange = sort(unique(S.xLGrid));
            
            S = surfaces{ii};
            subplot(1, length(surfaces), ii)
            hold on
                        
            % > Surface
            IResiduals = S.IResiduals;
            gridMin = min(IResiduals(:));
            gridMax = max(IResiduals(:));
            
            if isfield(P.data, 'stddevMSE')
                stddevMSE = P.data.stddevMSE;
            else
                stddevMSE = 20/100 * gridMin;  % Default to 20%.
            end
            isWithin1SD = (abs(IResiduals - gridMin) <= stddevMSE);
            isWithin3SD = ~isWithin1SD & (abs(IResiduals - gridMin) <= 3*stddevMSE);
            
            CO(:,:,1) = isWithin3SD * 0.5; %  red
            CO(:,:,2) = (isWithin1SD|isWithin3SD) * 0.5; % green
            CO(:,:,3) = ~(isWithin1SD|isWithin3SD) .* (gridMax-IResiduals)/gridMax; % blue
            caxis([gridMin, gridMin + stddevMSE]);
            
            surf(xLRange, nLRange, IResiduals, CO,...
                'HandleVisibility', 'off', ...
                'FaceColor', 'interp');
            
            % > Contour
            numLevels = 15;
            levels = logspace(log10(min(S.IResiduals(:))), log10(max(S.IResiduals(:))), numLevels); % non-linear spacing
%             levels = linspace(min(S.IResiduals(:)), max(S.IResiduals(:)), numLevels); % linear spacing
            contour3(xLRange, nLRange, S.IResiduals, ...
                levels, ...
                'Color', 'r', ...
                'HandleVisibility', 'off');
                        
%             dim = [0.2 0.5 0.3 0.3];
%             txt = sprintf("SD = %.2f, %.1f%% of minimum MSE", ...
%                 P.data.stddevMSE, P.data.stddevMSE/gridMin*100);            
%             annotation('textbox', dim, 'String', txt, ...
%                 'FitBoxToText', 'on', ...
%                 'BackgroundColor', 'white');
            
            title(sprintf("%s: %s", P.patientCode, S.name))
            
            xlabel("$x_L$ [1/min]")
            ylabel("$n_L$ [-]")
            zlabel("Mean of squared errors [(mU/min)^2]")
        end
        
    end
    
sprintf("%s \n%.2f \n\%d \n.1f%%\n\n", ...
    P.patientCode, P.data.stddevMSE, gridMin, P.data.stddevMSE/gridMin*100)
end

end

%% Functions
function [IResiduals, simI] = EvaluateGrid(PArray, nLGrid, xLGrid, savename)
global C

global FILEFORMAT

ISimulated = cell(size(nLGrid));
IResiduals = zeros(size(nLGrid));
for ii = 1:numel(nLGrid)
    if length(PArray) == 1
        copyP = PArray(1);
    else
        copyP = PArray(ii);
    end
    
    fprintf('\n%s: Trialling nL/xL = %g/%g in forward simulation (%d/%d). ', ...
        copyP.patientCode, nLGrid(ii), xLGrid(ii), ii, numel(nLGrid))
    
    % Apply nL/xL for iteration.
    copyP.results.nL = nLGrid(ii) * ones(size(copyP.results.tArray));
    copyP.results.xL = xLGrid(ii) * ones(size(copyP.results.tArray));
    
    % Get other parameters and forward simulate models.
    copyP = FindGutEmptyingRate(copyP);
    copyP = FitInsulinSensitivity(copyP, false);
    copyP = SolveSystem(copyP);
    
    % Determine error.
    if (copyP.source == "Detemir")
        [tI, vI] = GetSimTime(copyP, copyP.data.ITotal);        % Data [pmol/L]
        simI = C.mU2pmol(copyP.results.I + copyP.results.IDF);  % Sim [mU/L] -> [pmol/L]
    else
        [tI, vI] = GetSimTime(copyP, copyP.data.I);             % Data [pmol/L]
        simI = C.mU2pmol(copyP.results.I);                      % Sim [mU/L] -> [pmol/L]
    end
    iiI = GetTimeIndex(tI, copyP.results.tArray);
    simI = simI(iiI);
    IErrors = simI - vI;
    
    % Save residuals.
    IResiduals(ii) = sum(IErrors.^2)/numel(vI);  % Mean Squared Errors
    ISimulated{ii} = simI;
    
    EstimateTimeRemaining(ii, numel(nLGrid))
end

% Export results.
save(ResultsPath(sprintf(FILEFORMAT, PArray.patientCode, savename)), ...
    'nLGrid', 'xLGrid', 'IResiduals', 'ISimulated')

end
