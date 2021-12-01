function P = SplineSim(P, xL)
% Recipe for fitting parameters to a model with a fixed xL and splines for
% nL.
% INPUTS:
%   P  - patient struct
%   xL - (optional) choose xL to fix
% OUTPUT:
%   P  - updated patient struct

defaultxL = 0.6;

%% Functions
P = EstimateInsulinSecretion(P);

% Fix xL to hard-code value.
if ~exist("xL", "var")
    xL = defaultxL;
end
P.results.xL = xL;


% Fit nL with splines over range.
P = FitSplinesnL(P);


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
% PlotGlucosenL(P);

end

