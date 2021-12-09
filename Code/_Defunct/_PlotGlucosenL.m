function PlotGlucosenL(P)

persistent nLGlucosePlot;

DP = DebugPlots().FitSplines;
if DP.nLGlucose
    
    if isempty(nLGlucosePlot)
        P0 = P;
        P0.patientNum = 999;
        nLGlucosePlot = MakeDebugPlot("nL and G Correlation", P0, DP);
    else
        figure(nLGlucosePlot)
    end
    
    
    [tG, vG] = GetData(P.data.G); % [mmol/L]
    iiG = GetTimeIndex(tG, P.results.tArray);
    
    vGPlot = vG;
    
    nLPlot = P.results.nL(iiG);
    
%     % Scatter
%     sct = plot(vGNorm, nL, 'x');
%     sct.DisplayName = "P" + P.patientNum;
    
    % Linear Interpolation
    A(:,1) = vGPlot;
    A(:,2) = ones(size(vGPlot));
    b = nLPlot;
    theta = A\b;
    m = theta(1);
    c = theta(2);
    
    GArray = [min(vGPlot) max(vGPlot)];
    plt = plot(GArray, m*GArray+c);
    plt.DisplayName = P.patientCode;
    
    xlabel("Measured $G$ [mmol/L]")
    ylabel("Fit $n_L$ [1/min]")
    
    legend
    
end
end