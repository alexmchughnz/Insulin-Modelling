function PArray = InvestigateSplinesSim(P, xL)
% Recipe for fitting parameters to a model with a fixed xL and splines for
% nL.
% INPUTS:
%   P  - patient struct
%   xL - (optional) choose xL to fix
% OUTPUT:
%   P  - updated patient struct

defaultxL = 0.6;
PArray = {};

%% Setup
P = EstimateInsulinSecretion(P);

% Fix xL to hard-code value.
if ~exist("xL", "var")
    xL = defaultxL;
end
P.results.xL = xL;

% Find GFast.
[~, vG] = GetData(P.data.G);
P.data.GFast = min(vG);


%% Functions
% Iterate order with knots at data.
splineOptions.knotType = "location";
splineOptions.knots = P.data.I.time;

orderArray = 1:4;
for ii = 1:numel(orderArray)
    splineOptions.order = orderArray(ii);
    newP = TagPatientCode(P, "data | order = "+splineOptions.order);
    newP = FinishSimulation(newP, splineOptions);
    PArray = [PArray newP];
end


% Iterate number and order of fixed splines.
splineOptions.knotType = "amount";

numberArray = [2, 5, 10, 20];
orderArray = 1:4;
for nn = 1:numel(numberArray)
    splineOptions.knots = numberArray(nn);

    for kk = 1:numel(orderArray)
        splineOptions.order = orderArray(kk);
        newP = TagPatientCode(P, "number = " + splineOptions.knots + " | order = " + splineOptions.order);
        newP = FinishSimulation(newP, splineOptions);
        PArray = [PArray newP];
    end
end

end




function P = FinishSimulation(P, splineOptions)

% Fit nL with splines over range.
P = FitSplinesnL(P, splineOptions);


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
