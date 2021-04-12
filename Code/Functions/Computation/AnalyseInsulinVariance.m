function P = AnalyseInsulinVariance(P, stddev, N)
% Runs a Monte Carlo simulation, varying insulin data by a normal
% distributed scale factor with SD = stddev*data.
% INPUT:
%   P        - patient struct
%   stddev   - proportional
%   N        - number of trials to run
% OUTPUT:
%   P   - modified patient struct with nL and xL

DP = DebugPlots().AnalyseInsulinVariance;

%% Setup
tData = P.data.I.time;
MSE = zeros(1, N);

runtime = tic;

%% Simulate
for ii = 1:N
    % Randomly vary data according to normal distribution.
    randNums = Inf;
    while any(abs(randNums) > 3)  % Limit to within 3 SDs.
        randNums = randn(size(tData));
    end
    noiseFactors = 1 + stddev*randNums;
    
    % Vary correct data based on trial.
    copyP = P;
    
    if isfield(P, 'ITotal')
        [~, vITotal] = GetData(P.data.ITotal);
        trialITotal = noiseFactors .* vITotal;
        copyP.data.ITotal.value = trialITotal;
    else
        [~, vI] = GetData(P.data.I);
        trialI = noiseFactors .* vI;
        copyP.data.I.value = trialI;
    end
    
    % Forward simulate with varied data.
    MSE(ii) = GetSimError(copyP);
    
    % Print time.
    runtime = PrintTimeRemaining("AnalyseInsulinVariance", runtime, ii, N, P);
end

stddevError = std(MSE);
P.persistents.stddevMSE = stddevError;

message = sprintf("1 std. dev. of MSE is %g" ,stddevError);
PrintStatusUpdate(P, message, true);

%% ==================================================
%% Debug Plots
% Error
if DP.Error
    MakeDebugPlot("Insulin Error", P, DP);
    
    histogram(MSE, 50, ...
        'Normalization', 'probability');
    
    xlabel("Mean Squared Error")
    ylabel("Probability")
    
%     title(sprintf("P%d: Distribution of Model MSEs with Noise SD = %g*data (N = %d)", ...
%         P.patientNum, stddev, N))
    
    txt = sprintf("SD = %g", stddevError);
    xlimits = xlim;
    ylimits = ylim;
    text(0.9*xlimits(end), 0.9*ylimits(end), txt);
end

end


%% Functions
function MSE = GetSimError(P)
% Get other parameters and forward simulate models.
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, false);

% Determine error.
[tI, vI] = GetData(P.data.I);  % Data [mU/L]
[~, simI] = GetResultsSample(P, tI, P.results.I);

error = simI - vI;
MSE = mean(error.^2);
end


