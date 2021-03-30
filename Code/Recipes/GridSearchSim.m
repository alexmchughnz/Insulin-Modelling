function PArray = GridSearchSim(P)
% Recipe for basic fitting and forward simulating a patient.
% nL/xL found by grid search.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

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

integralP = FitHepaticClearance(P);  % To get A and b matrices.
integralP.patientCode = integralP.patientCode + "(integral)";
P.results.integrals = integralP.results.integrals;
    
if ~HasPersistent(P, "stddevMSE")    
    stddev = 5/100;
    N = 1000;
    
    P = AnalyseInsulinVariance(P, stddev, N);
end

gridSettings = {[0 0.4], [0.3 0.95], 0.02};
newGrid = true;
P = FindOptimalHepaticClearance(P, newGrid, gridSettings{:});

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

PArray = {P integralP};

end

