function plots = DebugPlots(newPlots)

persistent DEBUGPLOTS;

if exist("newPlots", "var")
    DEBUGPLOTS = newPlots;
    
elseif isempty(DEBUGPLOTS)
    DEBUGPLOTS.FIGURES = []; 
    
    %% Functions    
    DEBUGPLOTS.EstimateInsulinSecretion.Uen = true;
    DEBUGPLOTS.EstimateInsulinSecretion.CPep = false;
    
    DEBUGPLOTS.MakeSplineBasisFunctions.Splines = false;
    
    DEBUGPLOTS.IntegralFit.GraphicalID = false;
    DEBUGPLOTS.IntegralFit.Convergence = false;
    
    DEBUGPLOTS.FitSplines.Splines = true; 
    DEBUGPLOTS.FitSplines.nLGlucose = false; 
    
    DEBUGPLOTS.GridSearchParameters.ErrorSurface = true;
    
    DEBUGPLOTS.AnalyseInsulinVariance.Error = true;
    
    DEBUGPLOTS.SolveSystem.Glucose = true;
    DEBUGPLOTS.SolveSystem.PlasmaInsulin = true;
    DEBUGPLOTS.SolveSystem.InterstitialInsulin = false;
    DEBUGPLOTS.SolveSystem.GlucoseComponents = true;
    DEBUGPLOTS.SolveSystem.CoefficientShapes = false; 
    
    %% Recipes
    DEBUGPLOTS.IterateParametersSim.SSEHeatmap = true;
end

plots = DEBUGPLOTS;

end