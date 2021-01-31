function [] = PrintTimeRemaining(counter, max, P)

    fitTime = duration(seconds(toc));
    timeRemaining = datestr(fitTime*(max - counter + 1), 'HH:MM:SS');
    
    
    message = sprintf('(%d/%d) %s remaining.', ...
        counter, max, timeRemaining);
    PrintStatusUpdate(P, message, true);   
    
    tic;  
end

