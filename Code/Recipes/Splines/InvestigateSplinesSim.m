function PArray = InvestigateSplinesSim(P, xL)
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
splineOptions.knotType = "location";
splineOptions.knots = P.data.I.time;

orderArray = 1:4;
for ii = 1:numel(orderArray)
    splineOptions.order = orderArray(ii);
    newP = TagPatientCode(P, "data | order = "+splineOptions.order);
    newP = SCLossSplineSim(newP, splineOptions);
    PArray{end+1} = newP;
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
        newP = SCLossSplineSim(newP, splineOptions);
        PArray = [PArray newP];
    end
end

end


