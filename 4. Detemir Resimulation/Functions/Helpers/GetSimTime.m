function [time, value] = GetSimTime(P, datastruct)
    time = minutes(datastruct.time - P.simTime(1));
    
    inSimTime = (0 <= time) & (time < P.simDuration()); % [logical]
    
    time = time(inSimTime);    
    value = datastruct.value(inSimTime);    
end