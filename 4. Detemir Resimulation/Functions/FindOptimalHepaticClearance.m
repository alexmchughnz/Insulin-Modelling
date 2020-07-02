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
nLBest = nLGrid(iiOptimal);
xLBest = xLGrid(iiOptimal);

P.results.nL = nLBest * ones(size(tArray));
P.results.xL = xLBest * ones(size(tArray));


end
