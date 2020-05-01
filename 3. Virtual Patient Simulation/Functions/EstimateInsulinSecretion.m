function uEn = EstimateInsulinSecretion(G)
% Estimates pancreatic insulin secretion rate (Uen).
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen


global GC

%% Solving
if (G < GFast)
    uEn = GC.uMin;
elseif (G >= uMax)
    uEn = GC.uMax;
else
    uEn = GC.kSec*G + GC.kOff;
end

end

