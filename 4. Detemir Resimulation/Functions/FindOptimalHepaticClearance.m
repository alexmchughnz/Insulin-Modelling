function P = FindOptimalHepaticClearance(P, method)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'find' mode - very slow!
% INPUT:
%   P      - patient struct
%   method - 'find' to perform grid search
%            'load' to load previously-generated residuals data
% OUTPUT:
%   P   - modified patient struct with nL and xL

global C
global DEBUGPLOTS

%% Setup
tArray = P.results.tArray;
nLRange = 0 : 0.02 : 0.30;
xLRange = 0.3 : 0.1 : 1;
[nLGrid, xLGrid] = meshgrid(nLRange, xLRange);

if isequal(method, 'find')
    IResiduals = zeros(size(nLGrid));
    for ii = 1:numel(nLGrid)
        fprintf('\nP%d: Trialling nL/xL = %.2f/%.1f in forward simulation\n', ...
            P.patientNum, nLGrid(ii), xLGrid(ii))
        
        copyP = P;
        
        % Apply nL/xL for iteration.
        copyP.results.nL = nLGrid(ii) * ones(size(tArray));
        copyP.results.xL = xLGrid(ii) * ones(size(tArray));
        
        % Get other parameters and forward simulate models.
        copyP = FindGutEmptyingRate(copyP);
        copyP = FitInsulinSensitivity(copyP, true);
        copyP = SolveSystem(copyP);
        
        % Determine error.
        [tITotal, vITotal] = GetSimTime(copyP, copyP.data.ITotal);  % Data [pmol/L]
        iiITotal = GetTimeIndex(tITotal, tArray);
        
        simITotal = C.mU2pmol(copyP.results.I + copyP.results.IDF);  % Sim [mU/L] -> [pmol/L]
        simITotal = simITotal(iiITotal);
        
        ITotalError = 100*abs((simITotal - vITotal) ./ vITotal);
        
        % Save residuals.
        IResiduals(ii) = norm(ITotalError);
        save(sprintf('./Results/residualsP%d', P.patientNum), ...
            'nLGrid', 'xLGrid', 'IResiduals')
        
    end
elseif isequal(method, 'load')
    load(sprintf('./Results/residualsP%d', P.patientNum), ...
        'nLGrid', 'xLGrid', 'IResiduals');
end


%% Find Optimal nL/xL
iiOptimal = find(IResiduals == min(IResiduals(:)));
bestnL = nLGrid(iiOptimal);
bestxL = xLGrid(iiOptimal);

P.results.nL = bestnL * ones(size(tArray));
P.results.xL = bestxL * ones(size(tArray));


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
    
    plt = plot3(0, 1, min(IResiduals(:)), 'y*');
    plt.DisplayName = '$n_L/x_L = 0/1$';
    
    title(sprintf("P%d: Error Surface of (I+IDF) Fitting", P.patientNum))
    xlabel("$n_L$ [-]")
    ylabel("$x_L$ [1/min]")
    zlabel("2-norm of residuals, $\psi$ [mU/min]")
    
    legend
    
end

end
