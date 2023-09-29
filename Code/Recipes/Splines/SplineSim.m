function P = SplineSim(P, splineOptions, xL)
% Recipe for fitting parameters to a model with a fixed xL and splines for
% nL.
% INPUTS:
%   P  - patient struct
%   xL - (optional) choose xL to fix
% OUTPUT:
%   P  - updated patient struct

%% Setup

% Default splineOptions.
if ~exist("defaultSplineOptions", "var")
    defaultSplineOptions.knotType = "amount";
    defaultSplineOptions.knots = 20;
    defaultSplineOptions.order = 3;
    defaultSplineOptions.maxRate = 0.001;
    defaultSplineOptions.constrain = true;
    
    splineOptions = defaultSplineOptions;
end
P.results.splineOptions = splineOptions;

% Default xL.
if ~exist("xL", "var")
    defaultxL = 0.6;
    
    xL = defaultxL;
end
P.results.xL = xL;


%% Functions
P = EstimateInsulinSecretion(P);

% Fit nL with splines over range.
P = FitSplinesnL(P, splineOptions);


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


% Solve.
P = SolveSystem(P, true);

end

