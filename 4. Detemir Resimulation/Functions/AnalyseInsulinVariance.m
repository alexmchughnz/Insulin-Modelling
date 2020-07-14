function P = AnalyseInsulinVariance(P, stddev, N)
% Find optimal nL and xL, using grid search.
% Runs a LOT of forward simulations in 'find' mode - very slow!
% INPUT:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with nL and xL

global DEBUGPLOTS

%% Setup
nL = 0.22;
xL = 0.50;
P.results.nL = nL*ones(size(P.results.tArray));
P.results.xL = xL*ones(size(P.results.tArray));

[~, vITotal] = GetSimTime(P, P.data.ITotal);

%% Simulate
MSE = zeros(1, N);
scaleFactors = cell(1, N);

MSE(1) = GetSimError(P);
scaleFactors{1} = zeros(size(vITotal));
for ii = 2:N+1
    % Randomly vary data according to normal distribution.
    scaleFactors{ii} = 1 + stddev*randn(size(vITotal));
    trialITotal =  scaleFactors{ii} .* vITotal;
    
    copyP = P;
    copyP.data.ITotal.value = trialITotal;
    MSE(ii) = GetSimError(copyP);
    
    EstimateTimeRemaining(ii, N+1);
end

save(ResultsPath(sprintf("montecarlodata%gP%d", stddev, P.patientNum)))

%% Debug Plots
DP = DEBUGPLOTS.AnalyseInsulinVariance;

% Error
if DP.Error
    MakeDebugPlot(P, DP);    
    
    histogram(MSE, 50, ...
        'Normalization', 'probability');
    
    xlabel("Mean Squared Error")
    
    title(sprintf("Distribution of Model MSEs to Data + Normally-Distributed Noise with SD = %g%% (N = %d)", ...
        stddev*100, N))
end

stdError = std(MSE);
fprintf("1 std. dev. of MSE is %g\n", stdError)

end



function MSE = GetSimError(P)
global C

% Get other parameters and forward simulate models.
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P, true);
P = SolveSystem(P);

% Determine error.
[tITotal, vITotal] = GetSimTime(P, P.data.ITotal);
iiITotal = GetTimeIndex(tITotal, P.results.tArray);
simITotal = C.mU2pmol(P.results.I + P.results.IDF);  % Sim [mU/L] -> [pmol/L]
simITotal = simITotal(iiITotal);

error = simITotal - vITotal;
MSE = sum(error.^2)/length(vITotal);
end


