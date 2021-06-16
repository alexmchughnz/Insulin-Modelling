function P = SplineSim(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Plots
plots = DebugPlots();

DebugPlots(plots);


%% Functions
P = EstimateInsulinSecretion(P);

P.results.xL = 0.7
numSplines = 10;

% Manually solving here in the meantime.
P = FitSplinesJLKnL(P, numSplines);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

