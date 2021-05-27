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


[A, b] = AssembleIntegralSystemnLxL(P);
P.results.xL = 0.7;
P.results.nL = A(:,1) \ (b - A(:,2)*P.results.xL);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

