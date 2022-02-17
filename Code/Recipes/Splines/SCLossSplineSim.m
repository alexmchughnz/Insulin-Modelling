function P = SCLossSplineSim(P, splineOptions, xL)
% Recipe for fitting parameters to a model with a fixed xL and splines for nL.
% Identifies Lex, loss from subcutaneous injection.
% INPUTS:
%   P  - patient struct
%   xL - (optional) choose xL to fix
% OUTPUT:
%   P  - updated patient struct

%% Setup
if ~exist("splineOptions", "var")
    splineOptions = {}; 
end

% Default xL.
if ~exist("xL", "var")
    defaultxL = 0.6;
    
    xL = defaultxL;
end
P.results.xL = xL;


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

P = FitSplinesnL(P, splineOptions);


% Find optimal JLK.
JLKRange = 0.2 : 0.01 : 1.00;
P = LineSearchOptimum(P, "results.JLK", JLKRange, @InsulinError, @ApplyInsulinLossFactor);
P = ApplyInsulinLossFactor(P);


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

