function P = SimpleSim(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

                                                                                                                                             
%% Functions
P = EstimateInsulinSecretion(P);  % Fit Uen.

P = IntegralFitParameters(P, @AssembleIntegralSystemnLxL);  % Fit nL/xL, integral method.

P = FindGutEmptyingRate(P);  % Find d2.

GFastRange = [0.1 : 0.1 : 1] * P.data.GFast;
P = FindOptimalValue(P, "data.GFast", GFastRange, @GlucoseError, @FitInsulinSensitivity);

P = FitInsulinSensitivity(P);

P = SolveSystem(P, true);  % Simulate G, I, Q, etc.

end

