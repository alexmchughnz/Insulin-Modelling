function [time, value] = GetSimTime(P, datastruct)
    time = minutes(datastruct.time - P.simTime(1));
    value = datastruct.value;    
end