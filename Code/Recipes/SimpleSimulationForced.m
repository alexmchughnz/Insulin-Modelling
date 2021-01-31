function P = SimpleSimulationForced(P, forcenLxL)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
%   forcenLxL - [nL xL] value to set
% OUTPUT:
%   P  - updated patient struct

%% Plots
global DEBUGPLOTS
DEBUGPLOTS.EstimateInsulinSecretion.Uen = 1;
DEBUGPLOTS.FitInsulinSensitivity.SI = 1;
DEBUGPLOTS.SolveSystem.Glucose = 1;
DEBUGPLOTS.SolveSystem.Insulin = 1;

%% Functions
P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P, forcenLxL);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

