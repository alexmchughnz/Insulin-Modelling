function PArray = GridSearchSim(P)
% Recipe for basic fitting and forward simulating a patient.
% nL/xL found by grid search.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.EstimateInsulinSecretion.Uen = 0;
plots.EstimateInsulinSecretion.CPep = 0;
plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;

plots.FindOptimalHepaticClearance.ErrorSurface = 1;
plots.FindOptimalHepaticClearance.ErrorComparison = 1;

DebugPlots(plots);

%% Setup
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);  % To get A and b matrices.

if ~HasPersistent(P, "stddevMSE")
    % Currently irrelevant since error condition has changed.
    %     stddev = 5/100;
    %     N = 1000;
    %     P = AnalyseInsulinVariance(P, stddev, N);
    
    P.persistents.stddevMSE = 1000;
end

% Copy out P results from integral.
integralP = P;
integralP.patientCode = integralP.patientCode + "(integral)";

gridSettings = {[0 1], [0 1], 0.02};

newGrid = true;
P = FindOptimalHepaticClearance(P, newGrid, gridSettings{:});

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);


PArray = {P integralP};

end

