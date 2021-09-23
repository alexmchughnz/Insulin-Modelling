function P = FixedxLSplineSim(P)
% Recipe for
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Plots
plots = DebugPlots();

    plots.EstimateInsulinSecretion.Uen = false;
    plots.EstimateInsulinSecretion.CPep = false;
    
    plots.SolveSystem.CoefficientShapes = false; 
    
    plots.MakeSplineBasisFunctions.Splines = false;
    
    plots.FitSplines.nLGlucose = true;
    
DebugPlots(plots);


%% Functions
P = EstimateInsulinSecretion(P);

% Fix xL to hard-code value.
P.results.xL = 0.6;


% Fit nL with splines over range.
numKnots = numel(P.data.I.value) + 1;
P = FitSplinesnL(P, numKnots);


% Find optimal JLK.
JLKRange = 0.1 : 0.05 : 1.5;

    function P = ClearPForJLK(P) 
    P.results.P2 = zeros(size(P.results.tArray));
    P.results.d2 = 0;
    P.results.SI = 0;
    end

P = LineSearchOptimum(P, "results.JLK", JLKRange, @InsulinError, {@ClearPForJLK, @ApplyInsulinLossFactor});
P = ApplyInsulinLossFactor(P);  % Update


% Find GFast.
[~, vG] = GetData(P.data.G);
P.data.GFast = min(vG);


% Find d2.
lbHalfLife = 5;
ubHalfLife = 95;
halfLifeRange = 1 ./ linspace(1/ubHalfLife, 1/lbHalfLife, 20);
d2Range = log(2)./halfLifeRange;
P = LineSearchOptimum(P, "results.d2", d2Range, @GlucoseError, @FitInsulinSensitivity);


% Fit SI.
P = FitInsulinSensitivity(P);


% Find optimal ks3.
ks3Range = P.parameters.SC.ks3 * [0.1 : 0.2 : 1.5];
P = LineSearchOptimum(P, "parameters.SC.ks3", ks3Range, @GlucoseError, @AddTrialInputs);


% Fit SI.
P = FitInsulinSensitivity(P);


% Solve.
P = SolveSystem(P, true);
PlotGlucosenL(P);

end

