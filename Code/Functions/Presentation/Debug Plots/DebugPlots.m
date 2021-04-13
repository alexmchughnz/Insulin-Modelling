function plots = DebugPlots(newPlots)

persistent DEBUGPLOTS;

if exist("newPlots", "var")
    DEBUGPLOTS = newPlots;
    
elseif isempty(DEBUGPLOTS)    
    DEBUGPLOTS.EstimateInsulinSecretion.Uen = true;
    DEBUGPLOTS.EstimateInsulinSecretion.CPep = true;
    
    DEBUGPLOTS.FitHepaticClearance.GraphicalID = true; 
    DEBUGPLOTS.FitHepaticClearance.Convergence = true;
    
    DEBUGPLOTS.FindOptimalHepaticClearance.ErrorSurface = true;
    
    DEBUGPLOTS.AnalyseInsulinVariance.Error = true;
    
    DEBUGPLOTS.SolveSystem.Glucose = true;
    DEBUGPLOTS.SolveSystem.Insulin = true;
    
    DEBUGPLOTS.FIGURES = [];    
end

plots = DEBUGPLOTS;

end