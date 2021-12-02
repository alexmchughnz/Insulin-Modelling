function P = SplineSim(P)
% Recipe for fitting parameters to a model with a fixed xL and splines for
% nL.
% INPUTS:
%   P  - patient struct
%   xL - (optional) choose xL to fix
% OUTPUT:
%   P  - updated patient struct

defaultxL = 0.6;

% 1) Use measured C-Peptide to estimate endogenous secretion, Uen (van Cauter model).
P = EstimateInsulinSecretion(P);

% 2) Fix first pass hepatic clearance xL to hard-coded value.
P.results.xL = defaultxL;


% 3) Use splines to fit general hepatic clearance nL over time.
%    Optimal nL(t) minimises error of insulin prediction to measured data.
P = FitSplinesnL(P);


% 4) Just assume fasting glucose at start (reasonable for most trials).
[~, vG] = GetData(P.data.G);
P.data.GFast = vG(1);


% 5) Line search for optimum gut glucose absorption rate d2 that minimises glucose error.
lbHalfLife = 5;
ubHalfLife = 95;
halfLifeRange = 1 ./ linspace(1/ubHalfLife, 1/lbHalfLife, 20);
d2Range = log(2)./halfLifeRange;
P = LineSearchOptimum(P, "results.d2", d2Range, @GlucoseError, @FitInsulinSensitivity);


% 6) Use measured G and I, interpolated minute-to-minute, to estimate GI
P = FitInsulinSensitivity(P);


% 7) Now we have identified or selected all parameters needed to forward simulate the model.
%    This is all wrapped in this function:
P = SolveSystem(P, true);

end

