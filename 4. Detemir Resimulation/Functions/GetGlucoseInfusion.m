function P = GetGlucoseInfusion(P)
% Uses a patient's MANUAL data to determine additional glucose infusion.
% Values are hard-coded. This only applies to Patient 1.
% INPUTS:
%   P   - patient struct
%   t   - queried time
% OUTPUT:   
%   P   - modified patient struct with GInfusion

global C

t = 0 : P.simDuration();
P.GInfusion = zeros(size(t)); % By default, no infusion.

if (P.patientNum == 1)
    % Information about infusion.    
    MAGIC_DEXTROSE_NUMBER = 1.54;  % Assume this is some kind of "how much 
                                   % glucose from 5% dextrose" factor.
                                   
    duration = 12;          % Duration of infusion [hrs]
    duration = duration*60; % ''                   [min]
    
    startTime  = datetime('31/03/2017 05:15');
    preSimTime = minutes(P.simTime(1) - startTime); % How long infusion ran before sim [min]
    
    startTime = 0 - preSimTime;         % Start of infusion [min]
    endTime   = duration - preSimTime;  % End of infusion [min]
    
    % Return infusion data.
    iiInfusion = (startTime <= t) & (t < endTime); % 1 if infusion active [logical]
    P.GInfusion = iiInfusion * MAGIC_DEXTROSE_NUMBER/C.MGlucose;  % Glucose infusion over sim time [mol/min]
end

end
