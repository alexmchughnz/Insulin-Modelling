function P = nLxLGridSearchSim(P)
% Recipe for basic fitting and forward simulating a patient.
% nL/xL found by grid search.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Options
gridOptions.range = {[0 0.5], [0.4 1]};
gridOptions.step = [0.02];


%% Plots
plots = DebugPlots();
DebugPlots(plots);


%% Setup
P = EstimateInsulinSecretion(P);

% Find measure of variance due to insulin error for this patient.
[P, hasSSE] = GetPersistent(P, "stddevMSE");
if ~hasSSE
    stdDevPc = 5/100;
    N = 1000;
    P = AnalyseInsulinVariance(P, stdDevPc, N);
end

newGrid = true;
P = GridSearchParameters(P, @AssembleIntegralSystemnLxL, newGrid, gridOptions);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

