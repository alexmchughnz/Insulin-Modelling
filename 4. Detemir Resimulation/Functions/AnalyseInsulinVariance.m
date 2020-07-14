function P = AnalyseInsulinVariance(P, stddev, N)
% Runs a Monte Carlo simulation, varying insulin data by a normal
% distributed scale factor with SD = stddev*data.
% INPUT:
%   P        - patient struct
%   stddev   - proportional
%   N        - number of trials to run
% OUTPUT:
%   P   - modified patient struct with nL and xL

global DEBUGPLOTS

%% Setup
[~, vITotal] = GetSimTime(P, P.data.ITotal);
MSE = zeros(1, N);
scaleFactors = cell(1, N);

%% Simulate
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

save(ResultsPath(sprintf("montecarlodata%gx%d P%d", stddev, N, P.patientNum)))

% Save to patient.
stddevError = std(MSE);
P.results.stdError = stddevError;

%% ==================================================
%% Debug Plots
DP = DEBUGPLOTS.AnalyseInsulinVariance;

% Error
if DP.Error
    MakeDebugPlot(P, DP);    
    
    histogram(MSE, 50, ...
        'Normalization', 'probability');
    
    xlabel("Mean Squared Error")
    ylabel("Probability")
    
    title(sprintf("Distribution of Model MSEs to Data plus Normally-Distributed Noise with SD = %g%% (N = %d)", ...
        stddev*100, N))
end

fprintf("1 std. dev. of MSE is %g\n", stddevError)

end


%% Functions
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


