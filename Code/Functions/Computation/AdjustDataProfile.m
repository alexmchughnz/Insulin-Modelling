function P = AdjustDataProfile(P, pattern, varargin)
% desc
% INPUT:
%   P        - patient struct
%   stddev   - proportional
%   N        - number of trials to run
% OUTPUT:
%   P   - modified patient struct with nL and xL

DP = DebugPlots().AdjustDataProfile;

AUC = @(profile) trapz(P.results.tArray, profile);

%% Setup
path = varargin;
 dataProfile = getfield(P, path{:});

expandedPattern = repmat(pattern, 1, ceil(length(dataProfile)/length(pattern)));
expandedPattern = expandedPattern(1:length(dataProfile));
expandedPattern = expandedPattern(:);

newProfile = dataProfile .* expandedPattern;
P = setfield(P, path{:}, newProfile);


%% Optimise
originalAUC = AUC(dataProfile);
newAUC = AUC(newProfile);

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

