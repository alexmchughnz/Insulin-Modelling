function P = SimpleSim(P)
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

nLxL = [NaN, 0.7];
P = FixAndFitParameters(P, @AssembleIntegralSystemnLxL, nLxL);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

