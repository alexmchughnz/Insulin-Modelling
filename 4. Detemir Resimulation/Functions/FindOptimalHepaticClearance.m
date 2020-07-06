function P = FindOptimalHepaticClearance(P, method, varargin)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'find' mode - very slow!
% INPUT:
%   P        - patient struct
%   method   - 'find' to perform grid search
%              'load' to load previously-generated residuals data
%              'improve' to load data and iterate it to some precision
%   varargin - with 'load', the filename to load
%            - with 'improve', {1} the filename to load and improve
%                              {2} desired [nL, xL] grid precision
% OUTPUT:
%   P   - modified patient struct with nL and xL

load('config', 'RESULTPATH');

global DEBUGPLOTS
global FILEFORMAT
global resultsfile

RESIDUALSFORMAT = "residuals nL[%g %g]@%g xL[%g %g]@%g";
FILEFORMAT = '%sP%d.mat';
resultsfile = @(filename) fullfile(RESULTPATH, filename);

%% Setup
if isequal(method, 'find')
    nLDelta = 0.2;
    xLDelta = 0.2;
    
    nLBounds = [0 1];
    xLBounds = [-1 0];
    
    nLRange = nLBounds(1) : nLDelta : nLBounds(end);
    xLRange = xLBounds(1) : xLDelta : xLBounds(end);
    [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
    
    savename = sprintf(RESIDUALSFORMAT, ...
        nLBounds, nLDelta, xLBounds, xLDelta);
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename);
    
elseif isequal(method, 'load')
    % Load by name.
    loadname = varargin{1};
    load(resultsfile(sprintf(FILEFORMAT, loadname, P.patientNum)), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    nLRange = nLGrid(1, :);
    xLRange = xLGrid(:, 1);
    
    
elseif isequal(method, 'improve')
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
        
        % Shift window without further precision if best value on boundary.
        if (iinLOpt == 1) || (iinLOpt == numSegments) ...
                || (iixLOpt == 1) || (iixLOpt == numSegments)
            
            % Do nothing and continue!
            
        else
            if (nLDelta == nLPrecision) && (xLDelta == xLPrecision)
                % We're at desired precision, and optimum isn't on edge.
                % Time to finish up.
                loop = false;
            end
            
            % Update values for more precision.
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
iiOptimal = find(IResiduals == min(IResiduals(:)));
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


function IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename)
global C

global FILEFORMAT
global resultsfile

IResiduals = zeros(size(nLGrid));
for ii = 1:numel(nLGrid)
    fprintf('\nP%d: Trialling nL/xL = %g/%g in forward simulation (%d/%d)\n', ...
        P.patientNum, nLGrid(ii), xLGrid(ii), ii, numel(nLGrid))
    
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
end

end
