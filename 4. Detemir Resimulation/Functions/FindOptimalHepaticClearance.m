function P = FindOptimalHepaticClearance(P, method, arg)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'find' mode - very slow!
% INPUT:
%   P      - patient struct
%   method - 'find' to perform grid search
%            'load' to load previously-generated residuals data
%            'improve' to load data and iterate it to some precision
%   arg    - with 'load', the digit of improved residuals to load (eg. 0)
%            with 'improve', the desired [nL, xL] grid precision
% OUTPUT:
%   P   - modified patient struct with nL and xL

global DEBUGPLOTS

%% Setup
if isequal(method, 'find')
    nLDelta = 0.05;
    xLDelta = 0.1;
    
    nLBounds = [0 1];
    xLBounds = [0 1];
    
    nLRange = nLBounds(1) : nLDelta : nLBounds(end);
    xLRange = xLBounds(1) : xLDelta : xLBounds(end);
    [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
    
    savename = sprintf("residuals%g%g%g%g", nLBounds, xLBounds);
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename);
    
elseif isequal(method, 'load')
    if ~exist('arg', 'var')
        savename = arg;
        load(sprintf('./Results/%sP%d', savename, P.patientNum), ...
            'nLGrid', 'xLGrid', 'IResiduals');
    else
        dataNum = arg;
        load(sprintf('./Results/improveresiduals%dP%d',dataNum, P.patientNum), ...
            'nLGrid', 'xLGrid', 'IResiduals');
    end
    nLRange = nLGrid(1, :);
    xLRange = xLGrid(:, 1);
    
    
elseif isequal(method, 'improve')
    nLPrecision = arg(1);
    xLPrecision = arg(2);
    
    load(sprintf('./Results/residualsP%d', P.patientNum), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    
    nLRange = nLGrid(1, :);
    xLRange = xLGrid(:, 1);    
    nLDelta = diff(nLRange(1:2));
    xLDelta = diff(xLRange(1:2));
    
    loop = 0;
    factor = 10;
    numSegments = 10;
    while (nLDelta > nLPrecision) && (xLDelta > nLPrecision)
        % Find optima of current grid.
        iiOpt = find(IResiduals == min(IResiduals(:)));
        [iixLOpt, iinLOpt] = ind2sub(size(IResiduals), iiOpt);
        nLOpt = nLRange(iinLOpt);
        xLOpt = xLRange(iixLOpt);
        
        % Set up new search range around optima.
        nLRange = linspace(nLOpt - nLDelta, nLOpt + nLDelta, numSegments);
        xLRange = linspace(xLOpt - xLDelta, xLOpt + xLDelta, numSegments);
        [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
        
        % Evaluate over grid.
        savename = sprintf('improveresiduals%d', loop);
        IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename);
        
        % Update values.
        if nLDelta > factor*nLPrecision
            nLDelta = nLDelta / factor;
        else
            nLDelta = nLPrecision;
        end
        
        if xLDelta > factor*xLPrecision
            xLDelta = xLDelta / factor;
        else
            xLDelta = xLPrecision;
        end
        
        loop = loop + 1;
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
    save(sprintf(['./Results/%sP%d'], savename, P.patientNum), ...
        'nLGrid', 'xLGrid', 'IResiduals')
end

end
