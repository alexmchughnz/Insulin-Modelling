function PArray = IterateAmountSplineSim(P)
%% Plots
plots = DebugPlots();

DebugPlots(plots);


%% Functions
P = EstimateInsulinSecretion(P);

numSplinesArray = 3:10;

for ii = 1:length(numSplinesArray)
    numSplines = numSplinesArray(ii);
    
    copyP = P;
    copyP.patientSuffix = sprintf("#splines=%d", numSplines);
    
    copyP = FitSplinesJLKxLnL(copyP, numSplines);
    
    copyP = FindGutEmptyingRate(copyP);
    copyP = FitInsulinSensitivity(copyP);
    copyP = SolveSystem(copyP, true);
    
    PArray{ii} = copyP;
end

end

