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
[tData, ~] = GetIFromITotal(P);
MSE = zeros(1, N);
scaleFactors = cell(1, N);

%% Simulate
for ii = 1:N
    % Randomly vary data according to normal distribution.
    randNums = randn(size(tData));
    while any(abs(randNums) > 3)  % Limit to within 3 SDs.
        randNums = randn(size(tData));
    end
    scaleFactors{ii} = 1 + stddev*randNums;
    
    % Vary correct data based on trial.
    copyP = P;
    
    if isfield(P, 'ITotal')
        [~, vITotal] = GetSimTime(P, P.data.ITotal);
        trialITotal = scaleFactors{ii} .* vITotal;
        copyP.data.ITotal.value = trialITotal;
    else
        [~, vI] = GetSimTime(P, P.data.I);
        trialI = scaleFactors{ii} .* vI;
        copyP.data.I.value = trialI;
    end
    
    % Forward simulate with varied data.
    MSE(ii) = GetSimError(copyP);
    
    EstimateTimeRemaining(ii, N);
end
stddevError = std(MSE);
save(ResultsPath(sprintf("%s_montecarlo%gx%d.mat", P.patientCode, stddev, N)))

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
    
    title(sprintf("P%d: Distribution of Model MSEs with Noise SD = %g*data (N = %d)", ...
        P.patientNum, stddev, N))
    
    txt = sprintf("SD = %g", stddevError);
    xlimits = xlim;
    ylimits = ylim;
    text(0.9*xlimits(end), 0.9*ylimits(end), txt);
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
if isfield(P, 'ITotal')
    [tI, vI] = GetSimTime(P, P.data.ITotal);  % Data [mU/L]
    simI = P.results.I + P.results.IDF;       % Sim [mU/L]
else
    [tI, vI] = GetSimTime(P, P.data.I);  % Data [mU/L]
    simI = P.results.I;                  % Sim [mU/L]
end

iiI = GetTimeIndex(tI, P.results.tArray);
simI = simI(iiI);

error = simI - vI;
MSE = sum(error.^2)/length(vI);
end


