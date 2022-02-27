function PArray = InvestigateConstraintsSim(P, xL)
% Recipe for fitting parameters to a model with a fixed xL and splines for
% nL.
% INPUTS:
%   P  - patient struct
%   xL - (optional) choose xL to fix
% OUTPUT:
%   P  - updated patient struct

PArray = {};

%% Functions

% Analyse constant nL.
splineOptions = {};
splineOptions.knotType = "amount";
splineOptions.knots = 2;
splineOptions.maxRate = 0;
splineOptions.order = 1;

newP = TagPatientCode(P, "max rate = 0");
newP = SCLossSplineSim(newP, splineOptions);
PArray{end+1} = newP;


% Iterate max rate.
splineOptions = {};
maxRateArray = 0.001 * 10.^[0:2];

for ii = 1:numel(maxRateArray)
    splineOptions.maxRate = maxRateArray(ii);
    newP = TagPatientCode(P, "max rate = " + splineOptions.maxRate);
    newP = SCLossSplineSim(newP, splineOptions);
    PArray{end+1} = newP;
end

PArray = AnalysenLGlucose(PArray);

end



