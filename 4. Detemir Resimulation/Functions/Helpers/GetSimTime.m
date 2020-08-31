function [time, value] = GetSimTime(P, datastruct)

    if isdatetime(P.data.simTime(1))
        time = minutes(datastruct.time - P.data.simTime(1));
    else
        time = datastruct.time;
    end
    
    inSimTime = (P.data.simTime(1) <= time) & (time <= P.data.simTime(end)); % [logical]
    
    time = time(inSimTime);    
    value = datastruct.value(inSimTime);    
end