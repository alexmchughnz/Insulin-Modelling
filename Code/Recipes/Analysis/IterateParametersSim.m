function POpt = IterateParametersSim(P)
% INPUTS:
%   P  - patient struct
% OUTPUT:
%   PArray  - updated patient structs

%% Plots
plots = DebugPlots();

plots.EstimateInsulinSecretion.Uen = false;
plots.EstimateInsulinSecretion.CPep = false;
plots.IntegralFit.GraphicalID = false;
plots.IntegralFit.Convergence = false;
plots.SolveSystem.Insulin = false;
plots.SolveSystem.Glucose = false;
plots.SolveSystem.CoefficientShapes = false;

DebugPlots(plots);

%% Setup
PArray = {};

%% Functions
JLKArray = 0.5 : 0.1 : 1;  % "Justified Loss from injection in Knee"
ks3RateArray = 1 : 1 : 5;

% Find measure of variance due to insulin error for this patient.
[P, hasSSE] = GetPersistent(P, "stddevSSE");
if ~hasSSE
    stdDevPc = 5/100;
    N = 1000;
    P = AnalyseInsulinVariance(P, stdDevPc, N);
end

for kk = 1:numel(ks3RateArray)
    % Apply change to ks3 parameter.
    kkP = P;
    kkP = ScalePatientField(kkP, ks3RateArray(kk), "parameters", "SC", "ks3");
    
    for jj = 1:numel(JLKArray)
        % Apply IInput proportion.
        jjP = ScalePatientField(kkP, JLKArray(jj), "data", "vIBolus");
        jjP = AddBolusArrays(jjP);
        
        forceReSim = true;
        jjP = AddPlasmaInsulinInputArray(jjP, forceReSim); % Re-simulate SC model to change QLocal.
        
        % Now run full simulation for patient.
        jjP = SimpleSim(jjP);
        
        % Save.
        PArray{end+1} = jjP;
    end
end

% Hacky to only display *optimal* insulin plot.
plots.SolveSystem.Insulin = true;
plots.SolveSystem.CoefficientShapes = true;
DebugPlots(plots);
[~, iiOpt] = min(cellfun(@(P) P.results.fits.insulinSSE, PArray));

POpt = PArray{iiOpt};
SolveSystem(POpt, true);
POpt.patientCode = P.patientCode;

%% Plotting
plotvars.P = P;
plotvars.JLKArray = JLKArray;
plotvars.ks3RateArray = ks3RateArray;
MakePlots(PArray, plotvars);

end


function MakePlots(PArray, plotvars)
DP = DebugPlots().IterateParametersSim;

%% SSE Heatmap
if DP.SSEHeatmap   
    % Get grid of SSE values.
    gridDim = [numel(plotvars.JLKArray), numel(plotvars.ks3RateArray)];
    
    SSEArray = cellfun(@(P) P.results.fits.insulinSSE, PArray);    
    SSEGrid = reshape(SSEArray, gridDim);
    
    % Remove non-optimal values from grid.
    stddevSSE = plotvars.P.persistents.stddevSSE;
    minSSE =  min(SSEArray(:));
    isOptimalZone = abs(SSEArray - minSSE) <= stddevSSE;
    
    optimalSSEGrid = SSEGrid;
    optimalSSEGrid(~isOptimalZone) = NaN;
    
    
    F = MakeDebugPlot("SSE Heatmap", plotvars.P, DP);
    hold off  % Needed for heatmap.
    
    plotTitle = F.Children(1).Title.String;
    
    xvalues = plotvars.ks3RateArray; 
    yvalues = plotvars.JLKArray;
    h = heatmap(xvalues, yvalues, optimalSSEGrid);
    
    h.Title = plotTitle;
    h.XLabel = "ks3 Multiplier";
    h.YLabel = "JLK (proportion of SC injection present)";
    
end
end


