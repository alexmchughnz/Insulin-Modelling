function P = SimpleSimulation(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

P = EstimateInsulinSecretion(P);
P = FitHepaticClearance(P);
P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P);

end

