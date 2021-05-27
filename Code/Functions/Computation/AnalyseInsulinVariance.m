function P = AnalyseInsulinVariance(P, stdDevPc, N)
% Runs a Monte Carlo simulation, varying insulin data by a normal
% distributed scale factor with SD = stddev*data.
% INPUT:
%   P        - patient struct
%   stddevPc - proportional
%   N        - number of trials to run
% OUTPUT:
%   P   - modified patient struct with nL and xL

%% Setup
tData = P.data.I.time;
MSE = zeros(1, N);

runtime = tic;

%% Simulate 
% Ensure parameters are present.
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);

for ii = 1:N
    % Randomly vary data according to normal distribution.
    randNums = Inf;
    while any(abs(randNums) > 3)  % Limit to within 3 SDs.
        randNums = randn(size(tData));
    end
    varianceMultipliers = 1 + stdDevPc*randNums;
    
    % Vary correct data based on trial.
    copyP = ScalePatientField(P, varianceMultipliers, "data", "I", "value");
    
    % Forward simulate with varied data.
    [SSE(ii), MSE(ii)] = GetSimErrors(copyP);
    
    % Print time.
    runtime = PrintTimeRemaining("AnalyseInsulinVariance", runtime, ii, N, P, false, 50);
end

P.persistents.stddevMSE = std(MSE);
P.persistents.stddevSSE = std(SSE);

message1 = sprintf("1 std. dev. of SSE is %g", P.persistents.stddevSSE);
message2 = sprintf("1 std. dev. of MSE is %g", P.persistents.stddevMSE);
PrintStatusUpdate(P, message1, true);
PrintStatusUpdate(P, message2, true);


%% Plotting
plotvars.stddev = stdDevPc;
plotvars.N = N;
plotvars.MSE = MSE;
MakePlots(P, plotvars);


end


function [SSE, MSE] = GetSimErrors(P)
% Forward simulate models.
P = SolveSystem(P, false);

% Determine error.
[tI, vI] = GetData(P.data.I);  % Data [mU/L]
[~, simI] = GetResultsSample(P, tI, P.results.I);

error = simI - vI;
SSE = sum(error.^2);
MSE = mean(error.^2);
end


function MakePlots(P, plotvars)
DP = DebugPlots().AnalyseInsulinVariance;

%% Error
if DP.Error
    plotTitle = sprintf("Distribution of Model MSEs with Noise SD = %g*data (N = %d)", plotvars.stddev, plotvars.N);
    MakeDebugPlot(plotTitle, P, DP);
    
    histogram(plotvars.MSE, 50, ...
        'Normalization', 'probability');
    
    xlabel("Mean Squared Error")
    ylabel("Probability")  
    
    txt = sprintf("SD = %g", std(plotvars.MSE));
    xlimits = xlim;
    ylimits = ylim;
    text(0.9*xlimits(end), 0.9*ylimits(end), txt);
end

end

