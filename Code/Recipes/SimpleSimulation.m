function P = SimpleSimulation(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Plots
global DEBUGPLOTS
DEBUGPLOTS.EstimateInsulinSecretion.Uen = 1;
DEBUGPLOTS.FitInsulinSensitivity.SI = 1;
DEBUGPLOTS.PlotModelResults.Glucose = 1;
DEBUGPLOTS.PlotModelResults.Insulin = 1;

%% Functions
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P);

end

