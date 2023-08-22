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


% Analyse constant nL.
splineOptions = {};
splineOptions.knotType = "amount";
splineOptions.knots = 2;
splineOptions.maxRate = 0;
splineOptions.order = 1;

tag = "max rate = 0";
PArray{end+1} = RunSim(P, splineOptions, tag);


% Iterate order with knots at data.
measTimes = P.data.I.time;

splineOptions = {};
splineOptions.knotType = "location";
splineOptions.knots = measTimes;

orderArray = 1:4;
for ii = 1:numel(orderArray)
    splineOptions.order = orderArray(ii);
    tag = "data | order = "+splineOptions.order;
    
    PArray{end+1} = RunSim(P, splineOptions, tag);
end

% Data, double density.
midpoints = measTimes(1:end-1) + diff(measTimes)/2;

splineOptions = {};
splineOptions.knotType = "location";
splineOptions.knots = sort([measTimes(:); midpoints(:)]);

for ii = 1:numel(orderArray)
    splineOptions.order = orderArray(ii);
    tag = "data (x2) | order = "+splineOptions.order;
    
    PArray{end+1} = RunSim(P, splineOptions, tag);
end


% Iterate number and order of fixed splines.
splineOptions = {};
splineOptions.knotType = "amount";

numberArray = [2, 5, 10, 20, 50];
orderArray = 1:4;
for nn = 1:numel(numberArray)
    splineOptions.knots = numberArray(nn);

    for kk = 1:numel(orderArray)
        splineOptions.order = orderArray(kk);
        tag = "number = " + splineOptions.knots + " | order = " + splineOptions.order;
        
        PArray{end+1} = RunSim(P, splineOptions, tag);
    end
end

end



function newP = RunSim(P, splineOptions, tag)
    newP = TagPatientCode(P, tag);
    newP = SCLossSplineSim(newP, splineOptions);
end


