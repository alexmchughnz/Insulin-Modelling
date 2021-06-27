function PArray = IteratexLSplineSim(P)
%% Plots
plots = DebugPlots();

DebugPlots(plots);


%% Functions
P = EstimateInsulinSecretion(P);

numSplines = 10;
xLArray = 0.6 : 0.05: 0.8;

for ii = 1:length(xLArray)
    xL = xLArray(ii);
    
    copyP = P;
    copyP.patientSuffix = sprintf("xL=%.2f", xL);
    copyP.results.xL = xL;
    
    copyP = FitSplinesJLKnL(copyP, numSplines);
    
    copyP = FindGutEmptyingRate(copyP);
    copyP = FitInsulinSensitivity(copyP);
    copyP = SolveSystem(copyP, true);
    
    PArray{ii} = copyP;
end

end

