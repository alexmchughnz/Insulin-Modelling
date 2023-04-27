function PArray = InvestigateConstraintsSim(P, xL)
% Recipe for fitting parameters to a model with a fixed xL and splines for
% nL.
% INPUTS:
%   P  - patient struct
%   xL - (optional) choose xL to fix
% OUTPUT:
%   P  - updated patient struct

PArray = {};

defaultSplineOptions.knotType = "location";
defaultSplineOptions.knots = P.data.I.time;
defaultSplineOptions.order = 2;

%% Functions

% Analyse constant nL.
splineOptions = defaultSplineOptions;
splineOptions.maxRate = 0;

newP = TagPatientCode(P, "max rate = 0");
newP = SCLossSplineSim(newP, splineOptions);
PArray{end+1} = newP;


% Iterate max rate.
splineOptions = defaultSplineOptions;
maxRateArray = 0.001 * 10.^[-1:2];

for ii = 1:numel(maxRateArray)
    splineOptions.maxRate = maxRateArray(ii);
    newP = TagPatientCode(P, "max rate = " + splineOptions.maxRate);
    newP = SCLossSplineSim(newP, splineOptions);
    PArray{end+1} = newP;
end

PArray = AnalysenLGlucose(PArray);

end



