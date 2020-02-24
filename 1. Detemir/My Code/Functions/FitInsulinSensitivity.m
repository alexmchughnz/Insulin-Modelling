function P = FitInsulinSensitivity(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
%   GC  - glycaemic control parameter set
% OUTPUT:
%   P   - modified patient struct with SI

intervalTime = 360;  % Time per interval [min]
numIntervals = floor(P.trialLength/intervalTime);




P.SI = zeros(numIntervals, 1);



% Find indices of trial times in CGM data.
iiStart = find(P.G{3}.time == P.trialTime(1));   % Trial start [index]
iiEnd   = find(P.G{3}.time == P.trialTime(end)); % Trial end [index]

% Get value/time arrays over trial period.
tG = minutes(P.G{3}.time(iiStart:iiEnd) - P.G{3}.time(iiStart));
vG = P.G{3}.value(iiStart:iiEnd);
tI = minutes(P.I.time - P.I.time(1));
vI = P.I.value;

% Interpolate with a piecewise polynomial.
ppG = griddedInterpolant(tG, vG);
ppI = griddedInterpolant(tI, vI);


end


function [dydt] = GIModel(t, y, sys, ppG, ppI)
% Adapted from fitSI/FAERIES_integrals.
qsto1  = y(1);
qsto2  = y(2);
qgut   = y(3);
Q      = y(4);


dqsto1 = D(t) - GI.k21*qsto1;
dqsto2 = GI.k21*qsto1 - ;
dqgut


yd(1) = 

end
