function PArray = GridSearchLoad(P)
% Recipe for loading a previously-generate grid search.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

DebugPlots(plots);


%% Setup
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);  % To get A and b matrices.

% Copy out P results from integral.
integralP = P;
integralP.patientCode = integralP.patientCode + "(integral)";

newGrid = false;
P = FindOptimalHepaticClearance(P, newGrid);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);


PArray = {P integralP};

end

