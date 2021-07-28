function P = SplineSim(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Plots
plots = DebugPlots();

    plots.EstimateInsulinSecretion.Uen = false;
    plots.EstimateInsulinSecretion.CPep = false;
    
    plots.SolveSystem.CoefficientShapes = false; 
    
    plots.MakeSplineBasisFunctions.Splines = true;
    
    plots.FitSplines.nLGlucose = true;
    
DebugPlots(plots);


%% Functions
P = EstimateInsulinSecretion(P);

% Manually solving here in the meantime.
numKnots = numel(P.data.I.value) + 1;
P = FitSplinesJLKxLnL(P, numKnots);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

