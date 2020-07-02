function P = FindOptimalHepaticClearance(P)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations - very slow!
% INPUT:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with nL and xL

global C
global DEBUGPLOTS

%% Setup
nLRange = 0 : 0.02 : 0.30;
xLRange = 0.3 : 0.1 : 1;
[nLGrid, xLGrid] = meshgrid(nLRange, xLRange);

IResiduals = zeros(size(nLGrid));

tArray = P.results.tArray;

for ii = 1:numel(nLGrid)
    fprintf('P%d: Trialling nL/xL = %.2f/%.1f in forward simulation\n', ...
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
end

save(sprintf('residualsP%d', P.patientNum), 'nLGrid', 'xLGrid', 'IResiduals')


%% Find Optimal nL/xL
iiOptimal = find(IResiduals == min(IResiduals));
nLBest = nLGrid(iiOptimal);
xLBest = xLGrid(iiOptimal);

P.results.nL = nLBest * ones(size(tArray));
P.results.xL = xLBest * ones(size(tArray));


end
