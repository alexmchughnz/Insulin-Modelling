function P = AdjustDataProfile(P, patternTime, patternScale, varargin)
% desc
% INPUT:
%   P        - patient struct
%   stddev   - proportional
%   N        - number of trials to run
% OUTPUT:
%   P   - modified patient struct with nL and xL

DP = DebugPlots().AdjustDataProfile;

AUC = @(time, data) trapz(time, data);

%% Setup
path = varargin;
dataProfile = getfield(P, path{:});

%% Adjust
% Collect the Uen points at measurement times.  
ppData = griddedInterpolant(P.results.tArray, dataProfile);

measurementTimes = P.data.CPep.time;
deltaPatternTime = diff(patternTime(1:2)); % Change in t between desired points of adjustment.

% If new points are being added to profile, determine indices and
% interpolate.
if deltaPatternTime < 1
    iiMeasurement = GetTimeIndex(measurementTimes, P.results.tArray);
    
    % Find the indices at each multiple of deltaPatternTime between each
    % pair of measurement times.
    midpoints = iiMeasurement(1:end-1) + round(diff(iiMeasurement).*deltaPatternTime);
    newIndices = sort([iiMeasurement; midpoints]);
    
    newTime = P.results.tArray(newIndices);
else
    newTime = measurementTimes;
end

newDataPoints = ppData(newTime);

% Now, apply pattern to new data points.
expandedPattern = repmat(patternScale(:), ...
    ceil(length(newDataPoints)/length(patternScale)), 1);
expandedPattern = expandedPattern(1:length(newDataPoints));

adjustedData = newDataPoints .* expandedPattern;

% Finally, interpolate new data profile on a minute-wise grid for later
% calculations.
ppNewData = griddedInterpolant(newTime, adjustedData);
newDataProfile = ppNewData(P.results.tArray);

%% Optimise
originalAUC = AUC(P.results.tArray, dataProfile);
newAUC = AUC(P.results.tArray, newDataProfile);

message1 = sprintf("Original AUC: %.2f", originalAUC);
message2 = sprintf("New AUC: %.2f", newAUC);
PrintStatusUpdate(P, message1, true);
PrintStatusUpdate(P, message2, true);


%% Save
P = setfield(P, path{:}, newDataProfile);


%% ==================================================


%% Debug Plots
% Before and After
if DP.BeforeAfter
    MakeDebugPlot(path{end} + "Before VS After Adjustment", P, DP);
    
    plt = plot(P.results.tArray, dataProfile, 'k:');
    plt.DisplayName = path{end};
    
    plt = plot(P.results.tArray, newDataProfile, 'r');
    plt.DisplayName = "Adjusted " + path{end};
    
    xlabel("Time [min]")
    ylabel(path{end})
    
    legend()
end

end

