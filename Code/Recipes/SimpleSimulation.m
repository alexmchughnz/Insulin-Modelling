function P = SimpleSimulation(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Plots
plots = DebugPlots();

plots.EstimateInsulinSecretion.Uen = 1;
plots.EstimateInsulinSecretion.CPep = 1;
plots.FitHepaticClearance.GraphicalID = 1;
plots.FitInsulinSensitivity.SI = 1;
plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;

DebugPlots(plots);

%% Functions
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

