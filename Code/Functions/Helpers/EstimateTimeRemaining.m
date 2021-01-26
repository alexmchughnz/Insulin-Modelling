function [] = EstimateTimeRemaining(counter, max)
%ESTIMATETIMEREMAINING Summary of this function goes here
%   Detailed explanation goes here
    fitTime = duration(seconds(toc));
    timeRemaining = datestr(fitTime*(max - counter + 1), 'HH:MM:SS');
    
    
    fprintf('(%d/%d) Estimated time remaining: %s\n', ...
        counter, max, timeRemaining)
    
    tic;  
end

