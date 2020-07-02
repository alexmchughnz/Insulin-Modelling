function P = FindOptimalHepaticClearance(P, method, precision)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'find' mode - very slow!
% INPUT:
%   P      - patient struct
%   method - 'find' to perform grid search
%            'load' to load previously-generated residuals data
% OUTPUT:
%   P   - modified patient struct with nL and xL

global DEBUGPLOTS

%% Setup
nLDelta = 0.01;
xLDelta = 0.02;

nLRange = 0 : nLDelta : 0.30;
xLRange = 0.3 : xLDelta : 1;
[nLGrid, xLGrid] = meshgrid(nLRange, xLRange);

if isequal(method, 'find')
    IResiduals = EvaluateGrid(P, nLGrid, xLGrid, 'residuals');
    
elseif isequal(method, 'load')
    load(sprintf('./Results/residualsP%d', P.patientNum), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    
elseif isequal(method, 'improve')
    load(sprintf('./Results/residualsP%d', P.patientNum), ...
        'nLGrid', 'xLGrid', 'IResiduals');
    
    nLPrecision = precision(1);
    xLPrecision = precision(2);
    
    loop = 0;
    while (nLDelta > nLPrecision) && (xLDelta > nLPrecision)
        % Find optima of current grid.
        iiOpt = find(IResiduals == min(IResiduals(:)));
        [iixLOpt, iinLOpt] = ind2sub(size(IResiduals), iiOpt);
        nLOpt = nLRange(iinLOpt);
        xLOpt = xLRange(iixLOpt);
        
        % Set up new search range around optima.
        nLRange = linspace(nLOpt - nLDelta/2, nLOpt + nLDelta/2, 10);
        xLRange = linspace(xLOpt - xLDelta/2, xLOpt + xLDelta/2, 10);
        [nLGrid, xLGrid] = meshgrid(nLRange, xLRange);
        
        % Evaluate over grid.
        savename = sprintf('improveresiduals%d', loop);
        IResiduals = EvaluateGrid(P, nLGrid, xLGrid, savename);
        
        % Update values.
        if nLDelta > 10*nLPrecision
            nLDelta = nLDelta / 10;
        else
            nLDelta = nLPrecision;
        end
        
        if xLDelta > 10*xLPrecision
            xLDelta = xLDelta / 10;
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
    
    levels = logspace(log10(min(IResiduals(:))), log10(max(IResiduals(:))), 10);
    contour3(nLRange, xLRange, IResiduals, ...
        levels, ...
        'Color', 'r', ...
        'HandleVisibility', 'off');
    
    plt = plot3(bestnL, bestxL, min(IResiduals(:)), 'r*');
    plt.DisplayName = 'Optimal Point';
    
    plt = plot3(0, 1, IResiduals(end, 1), 'g*');
    plt.DisplayName = '$n_L/x_L = 0/1$';
    
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
    fprintf('\nP%d: Trialling nL/xL = %.2f/%.1f in forward simulation (%d/%d)\n', ...
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
