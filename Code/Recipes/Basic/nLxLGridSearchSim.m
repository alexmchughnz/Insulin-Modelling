function P = nLxLGridSearchSim(P)
% Recipe for basic fitting and forward simulating a patient.
% nL/xL found by grid search.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Options
gridOptions.range = {[0 1], [0 1]};
gridOptions.step = [0.2];


%% Plots
plots = DebugPlots();
DebugPlots(plots);


%% Setup
P = EstimateInsulinSecretion(P);

[P, hasMSE] = GetPersistent(P, "stddevMSE");
if ~hasMSE
    % Currently irrelevant since error condition has changed.
    %     stddev = 5/100;
    %     N = 1000;
    %     P = AnalyseInsulinVariance(P, stddev, N);
    
    P.persistents.stddevMSE = 1000;
end

newGrid = true;
P = GridSearchParameters(P, @AssembleIntegralSystemnLxL, newGrid, gridOptions);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);
P = SolveSystem(P, true);

end

