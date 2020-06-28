function [time, value] = GetSimTime(P, datastruct)
    time = minutes(datastruct.time - P.data.simTime(1));
    
    inSimTime = (0 <= time) & (time < P.data.simDuration()); % [logical]
    
    time = time(inSimTime);    
    value = datastruct.value(inSimTime);    
end