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
    
    plots.FitSplines.nLGlucose = false;
    
DebugPlots(plots);


%% Functions
P = EstimateInsulinSecretion(P);

numKnots = numel(P.data.I.value) + 1;
P = FitSplinesxLnL(P, numKnots);

halfLifeRange = 5 : 10 : 95;
d2Range = log(2)./halfLifeRange;
P = FindOptimalValue(P, "results.d2", d2Range, @GlucoseError, @FitInsulinSensitivity);

GFastRange = [0.1 : 0.1 : 1] * P.data.GFast;
P = FindOptimalValue(P, "data.GFast", GFastRange, @GlucoseError, @FitInsulinSensitivity);

ks3Range = [0.1 : 0.1 : 2] * P.parameters.SC.ks3;
P = FindOptimalValue(P, "parameters.SC.ks3", ks3Range, @InsulinError, @FitInsulinSensitivity);

P = FitSplinesJLKxLnL(P, numKnots);

P = FitInsulinSensitivity(P);

P = SolveSystem(P, true);

end

