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

load('config', 'RESULTPATH');

global DEBUGPLOTS
global FILEFORMAT
global resultsfile

GRIDFORMAT = "grid nL[%g %g]@%g xL[%g %g]@%g";
LINEFORMAT = "line nL=%g@%g to xL=%g@%g, t=%g";
LINE2DFORMAT = "2dline nL=%g@%g to xL=%g";
FILEFORMAT = '%sP%d.mat';
resultsfile = @(filename) fullfile(RESULTPATH, filename);

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
    
    [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
    
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
    
    [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
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
    load(resultsfile(sprintf(FILEFORMAT, loadname, P.patientNum)), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    
    words = split(loadname, ' ');
    type = words{1};
    if type == "line"
        % Reshape residuals onto grid.
        nLLine = nLGrid;
        nLRange = sort(unique(nLLine));
        xLLine = xLGrid;
        xLRange = sort(unique(xLLine));
        [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
        
        LineResiduals = IResiduals;
        
        IResiduals = nan(numel(xLRange), numel(nLRange));
        for ii = 1 : length(IResiduals)
            iinL = find(nLRange == nLLine(ii));
            iixL = find(xLRange == xLLine(ii));
            
            IResiduals(iixL, iinL) = LineResiduals(ii);
        end
        
    elseif type == "2dline"
        nLIntercept = nLGrid(xLGrid == 0);
        xLIntercept = xLGrid(nLGrid == 0);
        nLDelta = diff(nLGrid(1:2));
        
        nLtoxL = nLtoxLLineFun(nLIntercept, xLIntercept, nLDelta);
        
    else
        nLRange = nLGrid(1, :);
        xLRange = xLGrid(:, 1);
    end
    
    
elseif isequal(method, 'improve')
    type = 'grid';
    
    % Load by name.
    loadname = varargin{1};
    nLPrecision = varargin{2}(1);
    xLPrecision = varargin{2}(2);
    
    load(resultsfile(sprintf(FILEFORMAT, loadname, P.patientNum)), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    
    nLRange = nLGrid(1, :);
    xLRange = xLGrid(:, 1);
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
        [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
        
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
    
    load(resultsfile(sprintf(FILEFORMAT, loadname, P.patientNum)), ...
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
    bestIResiduals = EvaluateGrid(bestPArray, nLGrid, xLGrid, savename);
    
    savename = sprintf("%s worst", loadname);  % Worst Case
    worstIResiduals = EvaluateGrid(worstPArray, nLGrid, xLGrid, savename);
    
end

%% Find Optimal nL/xL
[minIResidual, iiOptimal] = min(IResiduals(:));
bestnL = nLGrid(iiOptimal);
bestxL = xLGrid(iiOptimal);

P.results.nL = bestnL * ones(size(P.results.tArray));
P.results.xL = bestxL * ones(size(P.results.tArray));

%% ------------------------------------------------------------------------

%% Debug Plots
DP = DEBUGPLOTS.FindOptimalHepaticClearance;

% Error Surface
if DP.ErrorSurface
    MakeDebugPlot(P, DP);
    hold on
    if type == "2dline"
        originalResiduals = IResiduals;
        plt = plot(nLGrid, originalResiduals);
        plt.DisplayName = 'Residuals';
        
        bestFile = sprintf(FILEFORMAT, loadname + " best", P.patientNum);
        if exist(bestFile, 'file')
            load(resultsfile(bestFile), ...
                'nLGrid', 'IResiduals');
            bestResiduals = IResiduals;
            plt = plot(nLGrid, bestResiduals);
            plt.DisplayName = 'Best Case Data Residuals';
        end
        
        worstFile = sprintf(FILEFORMAT, loadname + " worst", P.patientNum);
        if exist(worstFile, 'file')
            load(resultsfile(worstFile), ...
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
        
        figTitle = sprintf("P%d: Error Line of (I+IDF) Fitting Along $x_L = %.2f - %.2f n_L$", ...
            P.patientNum, xLIntercept, xLIntercept/nLIntercept);
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
            'name', 'Residuals');
        surfaces = [surfaces S];
        
        % Queue best/worst if they exist.
        bestFile = sprintf(FILEFORMAT, loadname + " best", P.patientNum);
        if exist(bestFile, 'file')
            load(resultsfile(bestFile), ...
                'nLGrid', 'IResiduals');
            S = struct('IResiduals', IResiduals, ...
                'nLGrid', nLGrid, ...
                'name', 'Best Case Data Residuals');
            surfaces = [surfaces S];
        end
        
        worstFile = sprintf(FILEFORMAT, loadname + " worst", P.patientNum);
        if exist(worstFile, 'file')
            load(resultsfile(worstFile), ...
                'nLGrid', 'IResiduals');
            S = struct('IResiduals', IResiduals, ...
                'nLGrid', nLGrid, ...
                'name', 'Worst Case Data Residuals');
            surfaces = [surfaces S];
        end
        
        % Plot each surface.
        for ii = 1:length(surfaces)
            subplot(1, length(surfaces), ii)
            hold on
            S = surfaces{ii};
            
            surf(nLRange, xLRange, S.IResiduals, ...
                'HandleVisibility', 'off');
            
            numLevels = 50;
            levels = logspace(log10(min(S.IResiduals(:))), log10(max(S.IResiduals(:))), numLevels);
            contour3(nLRange, xLRange, S.IResiduals, ...
                levels, ...
                'Color', 'r', ...
                'HandleVisibility', 'off');            
            
            
%             plt = plot3(bestnL, bestxL, min(IResiduals(:)), 'r*');
%             plt.DisplayName = 'Optimal Point';
%             
%             if ismember(0, nLRange) && ismember(1, xLRange)
%                 plt = plot3(0, 1, IResiduals(end, 1), 'g*');
%                 plt.DisplayName = '$n_L/x_L = 0/1$';
%             end

        title(sprintf("P%d: %s", P.patientNum, S.name))
        
        xlabel("$n_L$ [-]")
        ylabel("$x_L$ [1/min]")
        zlabel("2-norm of residuals, $\psi$ [mU/min]")
        end            
        
        tolerance = 2/100;
        bestFlatDomain = abs(surfaces{2}.IResiduals - IResiduals)./IResiduals <= tolerance;
        worstFlatDomain = abs(surfaces{3}.IResiduals - IResiduals)./IResiduals <= tolerance;
        flatDomain = bestFlatDomain & worstFlatDomain;
            subplot(1, length(surfaces), 1)
        plt = plot3(nLGrid, xLGrid, flatDomain*1000, 'y');
%             plt.DisplayName = 'Optimal Point';
        
    end
    
end

end

%% Functions
function [IResiduals, simITotal] = EvaluateGrid(PArray, nLGrid, xLGrid, savename)
global C

global FILEFORMAT
global resultsfile

ISimulated = cell(size(nLGrid));
IResiduals = zeros(size(nLGrid));
fitTime = duration();
for ii = 1:numel(nLGrid)
    tic
    fprintf('\nP%d: Trialling nL/xL = %g/%g in forward simulation (%d/%d). ', ...
        PArray.patientNum, nLGrid(ii), xLGrid(ii), ii, numel(nLGrid))
    fprintf('Estimated time remaining: %s\n', ...
        datestr(fitTime*(numel(nLGrid) - ii + 1),'HH:MM:SS'))
    
    if length(PArray) == 1
        copyP = PArray(1);
    else
        copyP = PArray(ii);
    end
    
    % Apply nL/xL for iteration.
    copyP.results.nL = nLGrid(ii) * ones(size(PArray.results.tArray));
    copyP.results.xL = xLGrid(ii) * ones(size(PArray.results.tArray));
    
    % Get other parameters and forward simulate models.
    copyP = FindGutEmptyingRate(copyP);
    copyP = FitInsulinSensitivity(copyP, true);
    copyP = SolveSystem(copyP);
    
    % Determine error.
    [tITotal, vITotal] = GetSimTime(copyP, copyP.data.ITotal);  % Data [pmol/L]
    iiITotal = GetTimeIndex(tITotal, PArray.results.tArray);
    
    simITotal = C.mU2pmol(copyP.results.I + copyP.results.IDF);  % Sim [mU/L] -> [pmol/L]
    simITotal = simITotal(iiITotal);
    
    ITotalError = 100*abs((simITotal - vITotal) ./ vITotal);
    
    % Save residuals.
    IResiduals(ii) = norm(ITotalError);
    ISimulated{ii} = simITotal;
    
    fitTime = duration(seconds(toc));
end

% Export results.
save(resultsfile(sprintf(FILEFORMAT, savename, PArray.patientNum)), ...
    'nLGrid', 'xLGrid', 'IResiduals', 'ISimulated')

end
