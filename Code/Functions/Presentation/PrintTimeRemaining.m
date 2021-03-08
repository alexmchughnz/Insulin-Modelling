function runtime = PrintTimeRemaining(activity, runtime, counter, max, P, showLine)
fitTime = duration(seconds(toc(runtime)));
timeRemaining = datestr(fitTime*(max - counter + 1), 'HH:MM:SS');

message = sprintf('(%d/%d) %s remaining for %s.', ...
    counter, max, timeRemaining, activity);

% Display

disp('')

if exist('showLine', 'var')
    if(showLine)
        PrintLine();
    end
end

PrintStatusUpdate(P, message, true);

if exist('showLine', 'var')
    if(showLine)
        PrintLine();
    end
end

runtime = tic;
end

