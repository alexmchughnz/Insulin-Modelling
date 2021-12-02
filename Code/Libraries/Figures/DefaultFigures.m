function figureList = DefaultFigures()

    %% Functions    
    figureList.EstimateInsulinSecretion.Uen = true;
    figureList.EstimateInsulinSecretion.CPep = true;
    
    figureList.MakeSplineBasisFunctions.BasisSplines = false;
    
    figureList.LineSearchOptimum.ErrorFunction = false;
    
    figureList.IntegralFit.GraphicalID = false;
    figureList.IntegralFit.Convergence = false;
    
    figureList.FitSplines.Splines = true; 
    figureList.FitSplines.nLGlucose = false; 
    
    figureList.GridSearchParameters.ErrorSurface = false;
    
    figureList.AnalyseInsulinVariance.Error = false;
    
    figureList.SolveSystem.PlasmaGlucose = true;
    figureList.SolveSystem.PlasmaInsulin = true;
    figureList.SolveSystem.InterstitialInsulin = false;
    figureList.SolveSystem.GlucoseComponents = false;
    figureList.SolveSystem.CoefficientShapes = false;
    
    figureList.MakeOGTTLui.Input = false;
    
    
    %% Recipes
    figureList.IterateParametersSim.SSEHeatmap = true;

end