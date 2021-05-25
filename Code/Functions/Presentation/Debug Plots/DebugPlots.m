function plots = DebugPlots(newPlots)

persistent DEBUGPLOTS;

if exist("newPlots", "var")
    DEBUGPLOTS = newPlots;
    
elseif isempty(DEBUGPLOTS)
    DEBUGPLOTS.FIGURES = []; 
    
    %% Functions    
    DEBUGPLOTS.EstimateInsulinSecretion.Uen = true;
    DEBUGPLOTS.EstimateInsulinSecretion.CPep = true;
    
    DEBUGPLOTS.IntegralFit.GraphicalID = true; 
    DEBUGPLOTS.IntegralFit.Convergence = true;
    
    DEBUGPLOTS.GridSearchParameters.ErrorSurface = true;
    
    DEBUGPLOTS.AnalyseInsulinVariance.Error = true;
    
    DEBUGPLOTS.SolveSystem.Glucose = true;
    DEBUGPLOTS.SolveSystem.Insulin = true;
    DEBUGPLOTS.SolveSystem.CoefficientShapes = true; 
    
    %% Recipes
    DEBUGPLOTS.IterateParametersSim.SSEHeatmap = true;
end

plots = DEBUGPLOTS;

end