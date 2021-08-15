function P = LineSearchJLKSim(P)
% Recipe for basic fitting and forward simulating a patient.
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Plots
plots = DebugPlots();

    plots.EstimateInsulinSecretion.Uen = false;
    plots.EstimateInsulinSecretion.CPep = false;
    
    plots.SolveSystem.CoefficientShapes = false; 
    
    plots.MakeSplineBasisFunctions.Splines = true;
    
    plots.FitSplines.nLGlucose = true;
    
DebugPlots(plots);


%% Functions
P = EstimateInsulinSecretion(P);

% Fit nL/xL.
numKnots = numel(P.data.I.value) + 1;
P = FitSplinesxLnL(P, numKnots);

P = FindGutEmptyingRate(P);
P = FitInsulinSensitivity(P);

P.results.xL

% Find optimal JLK.
[tI, vI] = GetData(P.data.I); % Measured I (for error comparison)

N = 50;
JLKGrid = linspace(1/100, 1, N);
for ii = 1:N
    JLK = JLKGrid(ii);
        
    message = "Searching at JLK = " + string(JLK) + "...";
    PrintStatusUpdate(P, message);     
    
    % Retrieve JLK value to simulate.
    copyP = ApplyInsulinLossFactor(P, JLK);
    
    % Simulate I(t, JLK) with resulting SI.
    copyP = SolveSystem(copyP);  
    
    % Find average error G(t, d2) to measured data.
    [~, simI] = GetResultsSample(copyP, tI, copyP.results.I);
    
    insulinSSE = (simI - vI).^2;
    IErrorGrid(ii) = mean(insulinSSE); 
end
[~, iiBest] = min(IErrorGrid);
P = ApplyInsulinLossFactor(P, JLKGrid(iiBest));


% Now re-fit nL/xL.
P = FitSplinesxLnL(P, numKnots);

JLKBest = JLKGrid(iiBest)


P.results.xL
P = SolveSystem(P, true);

end

