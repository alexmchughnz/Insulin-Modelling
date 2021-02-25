function runtime = PrintTimeRemaining(activity, runtime, counter, max, P)
    fitTime = duration(seconds(toc(runtime)));
    timeRemaining = datestr(fitTime*(max - counter + 1), 'HH:MM:SS');    
    
    message = sprintf('(%d/%d) %s remaining for %s.', ...
        counter, max, timeRemaining, activity);
    PrintStatusUpdate(P, message, true);   
    
    runtime = tic;  
end

