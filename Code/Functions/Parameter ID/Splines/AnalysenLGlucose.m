function PArray = AnalysenLGlucose(PArray)
% Plots, and stores ratio of, plasma glucose to nL clearance.

%% Setup
% Set up figure, and attach to first in array.
tag = "AnalysenLGlucose";
PArray{1} = AddFigure(PArray{1}, tag, "nLGlucose");
  

%% Functions
for ii = 1:numel(PArray)
    P = PArray{ii};
    
    % Get nL and G data points.
    [tG, vG] = GetData(P.data.G);
    nArray = SearchArray(tG, P.results.tArray);
    
    nLArray = P.results.nL(nArray);
    GArray = vG;
    
    % Linear regression.
    linear = 1;
    coeffs = polyfit(GArray, nLArray, linear);
    
    nLFit = polyval(coeffs, GArray);
    SSR = sum((nLFit - nLArray).^2);
    SST = (numel(nLArray)-1) * var(nLArray);
    Rsq = 1 - SSR/SST;
    
    % Save results.
    P.results.nLGlucose.coeffs = coeffs;
    P.results.nLGlucose.Rsq = Rsq;
    
    % Plot relationship on graph.
    plotvars.nLArray = nLArray;
    plotvars.GArray = GArray;
    AddToPlot(P, plotvars);
    
    PArray{ii} = P;
end

end


function AddToPlot(P, plotvars)
    % Raw Points
    plt = plot(plotvars.GArray, plotvars.nLArray, 'x');
    plt.HandleVisibility = 'off';
    color = plt.Color;
    
    % Linear Fit
    coeffs = P.results.nLGlucose.coeffs;
    
    nLFit = polyval(coeffs, plotvars.GArray);
    plt = plot(plotvars.GArray, nLFit);
    plt.DisplayName = P.patientCode;
    plt.Color = color;
    
    legend()
    
    ylim([0 0.5])
    
    xlabel("Plasma Glucose [mol/L]")
    ylabel("nL [1/min]")
end



