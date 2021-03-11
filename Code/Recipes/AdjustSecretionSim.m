function PArray = AdjustSecretionSim(P)
% Recipe for adjusting Uen and counter balancing by changing inputs.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.AdjustDataProfile.BeforeAfter = 1;
plots.EstimateInsulinSecretion.CPep = 1;
plots.SolveSystem.Glucose = 1;
plots.SolveSystem.Insulin = 1;
plots.FindOptimalHepaticClearance.ErrorSurface = 1;

DebugPlots(plots);


%% Setup
PArray = {};

P = EstimateInsulinSecretion(P);
fieldPath = ["results", "Uen"];

% Define the Uen adjustment strategies to use.
patterns{1} = 1;
patterns{2} = 1 + [-0.15 +0.15];

%% Functions
for pp = 1:length(patterns)
    pattern = patterns{pp};
    
    newP = AdjustDataProfile(P, pattern, fieldPath{:});
    [fitP, findP] = FitAndFindClearance(newP);
    PArray = [PArray fitP findP];
end

end




function [fitP, findP] = FitAndFindClearance(P)
% First, fit nL/xL and simulate.
fitP = P;
fitP.patientCode = fitP.patientCode + "fit";

fitP = FitHepaticClearance(fitP);
fitP = FindGutEmptyingRate(fitP);
fitP = FitInsulinSensitivity(fitP);
fitP = SolveSystem(fitP, true);

% Second, perform grid search.
findP = P;
findP.patientCode = findP.patientCode + "find";

if ~HasPersistent(findP, "stddevMSE")
    findP = FitHepaticClearance(findP);    
    stddev = 5/100;
    N = 1000;    
    findP = AnalyseInsulinVariance(findP, stddev, N);
end
gridSettings = {[-0.1 0.775], [0.075 0.95], 0.05};
findP = FindOptimalHepaticClearance(findP, "grid", gridSettings{:});

findP = FindGutEmptyingRate(findP);
findP = FitInsulinSensitivity(findP);
findP = SolveSystem(findP, true);
end
