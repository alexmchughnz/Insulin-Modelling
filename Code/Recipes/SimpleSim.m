function P = SimpleSim(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct


%% Plots
plots = DebugPlots();                      
    plots.EstimateInsulinSecretion.Uen = true;
    plots.EstimateInsulinSecretion.CPep = true;
    
    plots.SolveSystem.Glucose = true;
    plots.SolveSystem.PlasmaInsulin = true;
DebugPlots(plots);

%% Functions
P = EstimateInsulinSecretion(P);  % Fit Uen.

P = IntegralFitParameters(P, @AssembleIntegralSystemnLxL);  % Fit nL/xL, integral method.

% Find d2.
lbHalfLife = 5;
ubHalfLife = 95;
halfLifeRange = 1 ./ linspace(1/ubHalfLife, 1/lbHalfLife, 20);
d2Range = log(2)./halfLifeRange;
P = LineSearchOptimum(P, "results.d2", d2Range, @GlucoseError, @FitInsulinSensitivity);

P = FitInsulinSensitivity(P);

P = SolveSystem(P, true);  % Simulate G, I, Q, etc.

end

