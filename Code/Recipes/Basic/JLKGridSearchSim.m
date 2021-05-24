function P = JLKGridSearchSim(P)
% Recipe for basic fitting and forward simulating a patient.
% nL/xL found by grid search.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Options
gridOptions.range = {[0 1], [0 1]};
gridOptions.step = [0.1, 0.1];

%% Plots
plots = DebugPlots();

plots.EstimateInsulinSecretion.Uen = false;
plots.EstimateInsulinSecretion.CPep = false;
    
plots.FitHepaticClearance.GraphicalID = true; 
plots.FitHepaticClearance.Convergence = false;
    
plots.FindOptimalHepaticClearance.ErrorSurface = true;
    
plots.AnalyseInsulinVariance.Error = false;
    
plots.SolveSystem.Glucose = false;
DebugPlots(plots);


%% Setup
P = EstimateInsulinSecretion(P);
P.results.xL = 0.7; % FIXED

[P, hasMSE] = GetPersistent(P, "stddevMSE");
if ~hasMSE
%     Currently irrelevant since error condition has changed.
        stddev = 5/100;
        N = 1000;
        P = AnalyseInsulinVariance(P, stddev, N);
    
%     P.persistents.stddevMSE = 1000;
end

newGrid = true;
P = GridSearchParameters(P, @AssembleIntegralSystemJLKnL, newGrid, gridOptions);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

