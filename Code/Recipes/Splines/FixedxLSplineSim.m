function P = FixedxLSplineSim(P)
% Recipe for
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   P  - updated patient struct

%% Plots
plots = DebugPlots();

    plots.EstimateInsulinSecretion.Uen = false;
    plots.EstimateInsulinSecretion.CPep = false;
    
    plots.SolveSystem.CoefficientShapes = false; 
    
    plots.MakeSplineBasisFunctions.Splines = false;
    
    plots.FitSplines.nLGlucose = true;
    
DebugPlots(plots);


%% Functions
P = EstimateInsulinSecretion(P);

% Fix xL to hard-code value.
P.results.xL = 0.7;

% Fit nL with splines over range.
numKnots = numel(P.data.I.value) + 1;
P = FitSplinesnL(P, numKnots);

% Find d2 and fit SI.
halfLifeRange = 5 : 10 : 95;
d2Range = log(2)./halfLifeRange;
P = FindOptimalValue(P, "results.d2", d2Range, @GlucoseError, @FitInsulinSensitivity);

P = FitInsulinSensitivity(P);

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
JLKBest = JLKGrid(iiBest)
P = ApplyInsulinLossFactor(P, JLKBest);


P = SolveSystem(P, true);

end

