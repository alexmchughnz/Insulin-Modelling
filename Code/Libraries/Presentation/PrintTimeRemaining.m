function runtime = PrintTimeRemaining(activity, runtime, counter, max, ...
                                      P, showLine, interval)

% Settings
wrapFunc = @() 0;
if exist("showLine", "var")
    if(showLine)
        wrapFunc = @() PrintLine();
    end
end

if ~exist("interval", "var")
    interval = 1;
end

% Calculate
fitTime = duration(seconds(toc(runtime)));
timeRemaining = datestr(fitTime*(max - counter + 1), 'HH:MM:SS');

message = sprintf('(%d/%d) %s remaining for %s.', ...
    counter, max, timeRemaining, activity);

% Display
if mod(counter, interval) == 0
    disp('')
    disp('')
    disp('')
    
    wrapFunc();
    
    PrintStatusUpdate(P, message, true);
    
    wrapFunc();
end
    
runtime = tic;


end

