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
% Iterate order with knots at data.
maxRateArray = 0.001 * 10.^[0:2];
maxRateArray = [0 maxRateArray];
for ii = 1:numel(maxRateArray)
    splineOptions.maxRate = maxRateArray(ii);
    newP = TagPatientCode(P, "max rate = " + splineOptions.maxRate);
    newP = SCLossSplineSim(newP, splineOptions);
    PArray{end+1} = newP;
end

PArray = AnalysenLGlucose(PArray);

end



