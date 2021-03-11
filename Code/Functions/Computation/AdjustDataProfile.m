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

ppData = griddedInterpolant(P.results.tArray, dataProfile);

deltaIndex = diff(patternTime(1:2));
CPepTimes = P.data.CPep.time;

if deltaIndex < 1
iiUenPoints = GetTimeIndex(CPepTimes, P.results.tArray);

midpoints = iiUenPoints(1:end-1) + round(diff(iiUenPoints).*deltaIndex);
newIndices = sort([iiUenPoints; midpoints]);

newTime = P.results.tArray(newIndices);
else
    newTime = CPepTimes;
end
    
newData = ppData(newTime);

expandedPattern = repmat(patternScale(:), ...
                          ceil(length(newData)/length(patternScale)), 1);
expandedPattern = expandedPattern(1:length(newData));

adjustedData = newData .* expandedPattern;

ppNewData = griddedInterpolant(newTime, adjustedData);
newProfile = ppNewData(P.results.tArray);

P = setfield(P, path{:}, newProfile);


%% Optimise
originalAUC = AUC(P.results.tArray, dataProfile);
newAUC = AUC(P.results.tArray, newProfile);

message1 = sprintf("Original AUC: %.2f", originalAUC);
message2 = sprintf("New AUC: %.2f", newAUC);
PrintStatusUpdate(P, message1, true);
PrintStatusUpdate(P, message2, true);


%% ==================================================


%% Debug Plots
% Before and After
if DP.BeforeAfter    
    MakeDebugPlot(path{end} + "Before VS After Adjustment", P, DP);
    
    plt = plot(P.results.tArray, dataProfile, 'k:');
    plt.DisplayName = path{end};
    
    plt = plot(P.results.tArray, newProfile, 'r');
    plt.DisplayName = "Adjusted " + path{end};
    
    xlabel("Time [min]")
    ylabel(path{end})
    
    legend()
end

end

