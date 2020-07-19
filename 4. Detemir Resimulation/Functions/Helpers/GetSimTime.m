function [time, value] = GetSimTime(P, datastruct)

    if isdatetime(P.data.simTime(1))
        time = minutes(datastruct.time - P.data.simTime(1));
    else
        time = datastruct.time - P.data.simTime(1);
    end
    
    inSimTime = (0 <= time) & (time < P.data.simDuration()); % [logical]
    
    time = time(inSimTime);    
    value = datastruct.value(inSimTime);    
end