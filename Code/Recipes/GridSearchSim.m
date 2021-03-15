function P = GridSearchSim(P)
% Recipe for basic fitting and forward simulating a patient.
% nL/xL found by grid search.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Plots
plots = DebugPlots();

plots.EstimateInsulinSecretion.Uen = 1;
plots.EstimateInsulinSecretion.CPep = 1;
plots.FitInsulinSensitivity.SI = 1;
plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;

plots.FindOptimalHepaticClearance.ErrorSurface = 1;

DebugPlots(plots);

%% Setup
P = EstimateInsulinSecretion(P);

if ~HasPersistent(P, "stddevMSE")
    P = FitHepaticClearance(P);
    
    stddev = 5/100;
    N = 1000;
    
    P = AnalyseInsulinVariance(P, stddev, N);
end

gridSettings = {[0 0.4], [0.3 0.95], 0.025};
newGrid = true;
P = FindOptimalHepaticClearance(P, newGrid, gridSettings{:});

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

