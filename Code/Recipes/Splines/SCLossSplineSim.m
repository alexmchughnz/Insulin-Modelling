function P = SCLossSplineSim(P, xL)
% Recipe for fitting parameters to a model with a fixed xL and splines for nL.
% Identifies Lex, loss from subcutaneous injection.
% INPUTS:
%   P  - patient struct
%   xL - (optional) choose xL to fix
% OUTPUT:
%   P  - updated patient struct

defaultxL = 0.6;

%% Functions
P = EstimateInsulinSecretion(P);

% Scale Uen.
[~, nMaxI] = max(P.data.I.value);
tMaxI = P.data.I.time(nMaxI);

iiZero = SearchArray(0, P.results.tArray);
iiMaxI = SearchArray(tMaxI, P.results.tArray);

firstPhase = iiZero:iiMaxI;

P.results.Uen(firstPhase) = P.results.Uen(firstPhase) * 1.15;


% Fix xL to hard-code value.
if ~exist("xL", "var")
    xL = defaultxL;
end
P.results.xL = xL;


% Fit nL with splines over range.
splineOptions.knotType = "location";
splineOptions.knots = P.data.I.time;
splineOptions.order = 3;

P = FitSplinesnL(P, splineOptions);


% Find optimal JLK.
JLKRange = 0.2 : 0.01 : 1.00;
P = LineSearchOptimum(P, "results.JLK", JLKRange, @InsulinError, @ApplyInsulinLossFactor);
% P.results.JLK = 1;
P = ApplyInsulinLossFactor(P);


% % Find optimal ks3.
% ks3Range = P.parameters.SC.ks3 * [0.1 : 0.1 : 1.5];
% P = LineSearchOptimum(P, "parameters.SC.ks3", ks3Range, @InsulinError, @AddTrialInputs);


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

