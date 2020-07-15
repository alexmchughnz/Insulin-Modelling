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
for ii = 1:N
    % Randomly vary data according to normal distribution.
    randNums = randn(size(vITotal));
    while any(randNums > 3)  % Limit to within 3 SDs.
        randNums = randn(size(vITotal));
    end
    scaleFactors{ii} = 1 + stddev*randn(size(vITotal));
     
    trialITotal = scaleFactors{ii} .* vITotal;
    
    copyP = P;
    copyP.data.ITotal.value = trialITotal;
    MSE(ii) = GetSimError(copyP);
    
    EstimateTimeRemaining(ii, N);
end

save(ResultsPath(sprintf("montecarlodata%gx%d P%d.mat", stddev, N, P.patientNum)))

%% ==================================================
%% Debug Plots
DP = DEBUGPLOTS.AnalyseInsulinVariance;
stddevError = std(MSE);

% Error
if DP.Error
    MakeDebugPlot(P, DP);    
    
    histogram(MSE, 50, ...
        'Normalization', 'probability');
    
    xlabel("Mean Squared Error")
    ylabel("Probability")
    
    title(sprintf("P%d: Distribution of Model MSEs with Noise SD = %g*data (N = %d)", ...
        P.patientNum, stddev, N))
    
    txt = sprintf("SD = %g", stddevError);
    text(0, 1, txt);
end

fprintf("P%d: 1 std. dev. of MSE is %g\n", P.patientNum, stddevError)

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


