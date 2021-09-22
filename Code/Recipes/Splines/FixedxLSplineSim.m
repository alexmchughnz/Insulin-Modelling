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
P.results.xL = 0.6;


% Fit nL with splines over range.
numKnots = numel(P.data.I.value) + 1;
P = FitSplinesnL(P, numKnots);

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
    
    % Simulate I(t, JLK) with junk SI.
    copyP.results.P2 = zeros(size(P.results.tArray));
    copyP.results.d2 = 0;
    copyP.results.SI = 0;
    copyP = GCModel(copyP);  
    
    % Find average error G(t, d2) to measured data.
    [~, simI] = GetResultsSample(copyP, tI, copyP.results.I);
    
    insulinSSE = (simI - vI).^2;
    IErrorGrid(ii) = mean(insulinSSE); 
end
[~, iiBest] = min(IErrorGrid);
JLKBest = JLKGrid(iiBest);
P = ApplyInsulinLossFactor(P, JLKBest);

% Find GFast.
% GFastRange = 1 : 0.25 : 6;
% P = FindOptimalValue(P, "data.GFast", GFastRange, @GlucoseError, @FitInsulinSensitivity);
[~, vG] = GetData(P.data.G);
P.data.GFast = min(vG);

% Find d2.
lbHalfLife = 5;
ubHalfLife = 95;
halfLifeRange = 1 ./ linspace(1/ubHalfLife, 1/lbHalfLife, 20);
d2Range = log(2)./halfLifeRange;
P = FindOptimalValue(P, "results.d2", d2Range, @GlucoseError, @FitInsulinSensitivity);


% Fit SI.
P = FitInsulinSensitivity(P);


% Find optimal ks3.
ks3Range = P.parameters.SC.ks3 * [0.1 : 0.2 : 1.5];
funcsToApply = {@AddPlasmaInsulinInputArray @FitInsulinSensitivity};
P = FindOptimalValue(P, "parameters.SC.ks3", ks3Range, @GlucoseError, funcsToApply);


% Fit SI.
P = FitInsulinSensitivity(P);


% Solve.
P = SolveSystem(P, true);
PlotGlucosenL(P);

end

