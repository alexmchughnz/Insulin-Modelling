function PArray = AdjustSecretionSim(P)
% Recipe for adjusting Uen and counter balancing by changing inputs.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.AdjustDataProfile.BeforeAfter = 1;
plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;
plots.FindOptimalHepaticClearance.ErrorSurface = 1;

DebugPlots(plots);


%% Setup
PArray = {};

P = EstimateInsulinSecretion(P);
fieldPath = ["results", "Uen"];

% Define the Uen adjustment strategies to use.
% Time = Time points to apply pattern at. Intervals of <1 will add
% interpolated points.
% Scale = Pattern of scales to apply to points. Will be continued over
% length of data array.
MakePattern = @(time, scale) struct('time', time, 'scale', scale);

patterns{1} = MakePattern([1 2], ...
                          [1 1]);
patterns{2} = MakePattern([1 2], ...
                           1 + [-0.15 +0.15]);
patterns{3} = MakePattern([1 1.5 ], ...
                           1 + [-0.15 +0.15]);
                       
patterns{4} = MakePattern([1 2], ...
                           1 + [+0.15 -0.15]);
patterns{5} = MakePattern([1 1.5 ], ...
                           1 + [+0.15 -0.15]);

%% Functions
for pp = 1:length(patterns)
    pattern = patterns{pp};
    
    % Adjust patient metadata.
    newP = P;
    newP.patientNum = newP.patientNum*10 + pp;
    newP.patientCode = newP.patientCode + "." +string(pp);

    % Edit data profile, then fit clearance values.
    adjP = AdjustDataProfile(newP, pattern.time, pattern.scale, fieldPath{:});
    adjUen = adjP.results.Uen;
    
    [fitP, findP] = FitAndFindClearance(newP, adjUen);
    
    % Append simulated structs to PArray.
    PArray = [PArray fitP findP];
end

end




function [fitP, findP] = FitAndFindClearance(P, adjUen)
%% Fitting
fitP = P;
fitP.patientCode = fitP.patientCode + "fit";

% Apply adjusted Uen.
adjP = fitP;
adjP.results.Uen = adjUen;

% Fit.
adjP = FitHepaticClearance(adjP);
fitP.results.nL = adjP.results.nL;
fitP.results.xL = adjP.results.xL;

% Simulate.
fitP = FindGutEmptyingRate(fitP);
fitP = FitInsulinSensitivity(fitP);
fitP = SolveSystem(fitP, true);

%% Finding
findP = P;
findP.patientCode = findP.patientCode + "find";

if ~HasPersistent(findP, "stddevMSE")
    findP = FitHepaticClearance(findP);    
    stddev = 5/100;
    N = 1000;
    findP = AnalyseInsulinVariance(findP, stddev, N);
end

% Run grid search.
gridSettings = {[-0.1 0.775], [0.075 0.95], 0.025};
findP = FindOptimalHepaticClearance(findP, "grid", gridSettings{:});

% Simulate.
findP = FindGutEmptyingRate(findP);
findP = FitInsulinSensitivity(findP, false);
findP = SolveSystem(findP, true);
end
