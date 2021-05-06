function PArray = GridSearchSim(P)
% Recipe for basic fitting and forward simulating a patient.
% nL/xL found by grid search.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Settings
gridSettings = {[0.02 1], [0.02 1], 0.02};

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
P = FitHepaticClearance(P);  % To get A and b matrices.

[P, hasMSE] = GetPersistent(P, "stddevMSE");
if ~hasMSE
%     Currently irrelevant since error condition has changed.
        stddev = 5/100;
        N = 1000;
        P = AnalyseInsulinVariance(P, stddev, N);
    
%     P.persistents.stddevMSE = 1000;
end

% Copy out P results from integral.
integralP = P;
integralP.patientCode = integralP.patientCode + "(integral)";

newGrid = true;
P = FindOptimalHepaticClearance(P, newGrid, gridSettings{:});

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);


PArray = {P integralP};

end

