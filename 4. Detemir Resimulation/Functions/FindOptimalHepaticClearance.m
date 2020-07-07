function P = FindOptimalHepaticClearance(P, method, varargin)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'find' mode - very slow!
% INPUT:
%   P        - patient struct
%   method   - 'find' to perform grid search
%              'line' to perform grid search on line
%              '2dline' to perform 2D search on line
%              'load' to load previously-generated residuals data
%              'improve' to load data and iterate it to some precision
%   varargin - with 'line', [nL-intercept xL-intercept] to search over
%            - with 'load', the filename to load
%            - with 'improve', {1} the filename to load and improve
%                              {2} desired [nL, xL] grid precision
% OUTPUT:
%   P   - modified patient struct with nL and xL

load('config', 'RESULTPATH');

global DEBUGPLOTS
global FILEFORMAT
global resultsfile

RESIDUALSFORMAT = "residuals nL[%g %g]@%g xL[%g %g]@%g";
LINEFORMAT = "line nL=%g@%g to xL=%g@%g, t=%g";
LINE2DFORMAT = "2dline nL=%g@%g to xL=%g";
FILEFORMAT = '%sP%d.mat';
resultsfile = @(filename) fullfile(RESULTPATH, filename);

nLtoxLLineFun = @(nLIntercept, xLIntercept, nLDelta) ...
    (@(nL) RoundToMultiple(xLIntercept - xLIntercept/nLIntercept * nL, ...
    nLDelta));

%% Setup
if isequal(method, 'find')
    type = 'grid';
    
    nLDelta = 0.2;
    nLBounds = [0 1];
    nLRange = nLBounds(1) : nLDelta : nLBounds(end);
    
    xLDelta = 0.2;
    xLBounds = [-1 0];
    xLRange = xLBounds(1) : xLDelta : xLBounds(end);
    
    [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
    
    savename = sprintf(RESIDUALSFORMAT, ...
        nLBounds, nLDelta, xLBounds, xLDelta);
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename);
    
elseif isequal(method, 'line')
    type = 'line';
    
    nLIntercept = varargin {1}(1);
    xLIntercept = varargin {1}(2);
    
    % Establish full grid over area.
    nLDelta = 0.02;
    nLBounds = [0 nLIntercept];
    nLRange = nLBounds(1) : nLDelta : nLBounds(end);
    
    xLDelta = 0.02;
    xLBounds = [0 xLIntercept];
    xLRange = xLBounds(1) : xLDelta : xLBounds(end);
    
    [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
    
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
    IResiduals = nan(size(nLGrid));
    for ii = 1 : length(LineResiduals)
        iinL = find(nLRange == nLLine(ii));
        iixL = find(xLRange == xLLine(ii));
        
        IResiduals(iixL, iinL) = LineResiduals(ii);
    end
    
    
elseif isequal(method, '2dline')
    type = '2dline';
    
    nLIntercept = varargin {1}(1);
    xLIntercept = varargin {1}(2);
    
    % Establish full grid over area.
    nLDelta = 0.01;
    nLBounds = [0 nLIntercept];
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
        savename = sprintf(RESIDUALSFORMAT, ...
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
end

%% Find Optimal nL/xL
[minIResidual, iiOptimal] = min(IResiduals(:));
bestnL = nLGrid(iiOptimal);
bestxL = xLGrid(iiOptimal);

P.results.nL = bestnL * ones(size(P.results.tArray));
P.results.xL = bestxL * ones(size(P.results.tArray));


%% Debug Plots
DP = DEBUGPLOTS.FindOptimalHepaticClearance;

% Error Surface
if DP.ErrorSurface
    MakeDebugPlot(P, DP);
    hold on
    if type == "2dline"
        plt = plot(nLGrid, IResiduals);
        plt.DisplayName = 'Residuals';
        
        plt = plot(nLGrid(iiOptimal), minIResidual, 'r*');
        plt.DisplayName = 'Optimal Point';
        
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
        ylabel("2-norm of residuals, $\psi$ [mU/min]")
        
        grid on
        ax1.Position = ax2.Position;
    else
        
        surf(nLRange, xLRange, IResiduals, ...
            'HandleVisibility', 'off');
        
        numLevels = 50;
        levels = logspace(log10(min(IResiduals(:))), log10(max(IResiduals(:))), numLevels);
        contour3(nLRange, xLRange, IResiduals, ...
            levels, ...
            'Color', 'r', ...
            'HandleVisibility', 'off');
        
        plt = plot3(bestnL, bestxL, min(IResiduals(:)), 'r*');
        plt.DisplayName = 'Optimal Point';
        
        if ismember(0, nLRange) && ismember(1, xLRange)
            plt = plot3(0, 1, IResiduals(end, 1), 'g*');
            plt.DisplayName = '$n_L/x_L = 0/1$';
        end
        
        title(sprintf("P%d: Error Surface of (I+IDF) Fitting", P.patientNum))
        xlabel("$n_L$ [-]")
        ylabel("$x_L$ [1/min]")
        zlabel("2-norm of residuals, $\psi$ [mU/min]")
        
        legend
    end
    
end

end


function IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename)
global C

global FILEFORMAT
global resultsfile

IResiduals = zeros(size(nLGrid));
fitTime = duration();
for ii = 1:numel(nLGrid)
    tic
    fprintf('\nP%d: Trialling nL/xL = %g/%g in forward simulation (%d/%d). ', ...
        P.patientNum, nLGrid(ii), xLGrid(ii), ii, numel(nLGrid))
    fprintf('Estimated time remaining: %s\n', ...
        datestr(fitTime*(numel(nLGrid) - ii + 1),'HH:MM:SS'))
    
    copyP = P;
    
    % Apply nL/xL for iteration.
    copyP.results.nL = nLGrid(ii) * ones(size(P.results.tArray));
    copyP.results.xL = xLGrid(ii) * ones(size(P.results.tArray));
    
    % Get other parameters and forward simulate models.
    copyP = FindGutEmptyingRate(copyP);
    copyP = FitInsulinSensitivity(copyP, true);
    copyP = SolveSystem(copyP);
    
    % Determine error.
    [tITotal, vITotal] = GetSimTime(copyP, copyP.data.ITotal);  % Data [pmol/L]
    iiITotal = GetTimeIndex(tITotal, P.results.tArray);
    
    simITotal = C.mU2pmol(copyP.results.I + copyP.results.IDF);  % Sim [mU/L] -> [pmol/L]
    simITotal = simITotal(iiITotal);
    
    ITotalError = 100*abs((simITotal - vITotal) ./ vITotal);
    
    % Save residuals.
    IResiduals(ii) = norm(ITotalError);
    save(resultsfile(sprintf(FILEFORMAT, savename, P.patientNum)), ...
        'nLGrid', 'xLGrid', 'IResiduals')
    
    fitTime = duration(seconds(toc));
end

end
