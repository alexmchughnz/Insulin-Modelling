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
    isFakeData = (tData <= 0);
    randNums = Inf;
    while any(abs(randNums) > 3)  % Limit to within 3 SDs.
        randNums = randn(size(tData));
    end
    randNums(isFakeData) = randNums(1); % Apply same variance to fake data.
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
P = FitInsulinSensitivity(P, false);
P = SolveSystem(P, false);

% Determine error.
if isfield(P, 'ITotal')
    [tI, vI] = GetData(P.data.ITotal);  % Data [mU/L]
    simI = P.results.I + P.results.IDF;       % Sim [mU/L]
else
    [tI, vI] = GetData(P.data.I);  % Data [mU/L]
    simI = P.results.I;                  % Sim [mU/L]
end

iiI = GetTimeIndex(tI, P.results.tArray);
simI = simI(iiI);

error = simI - vI;
error = error(tI >= 0);  % Only evaluate error at true points.
MSE = mean(error.^2);
end


