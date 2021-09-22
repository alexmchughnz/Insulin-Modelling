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
    
    vGNorm = vG ./ P.data.GFast;
    
    nL = P.results.nL(iiG);
    
    % Scatter
    sct = plot(vGNorm, nL, 'x');
    sct.DisplayName = "P" + P.patientNum;
    
    % Linear Interpolation
    A(:,1) = vGNorm;
    A(:,2) = ones(size(vGNorm));
    b = nL;
    theta = A\b;
    m = theta(1);
    c = theta(2);
    
    GArray = 0:3;
    plt = plot(GArray, m*GArray+c, 'Color', sct.Color);
    plt.HandleVisibility = 'off';
    
    xlabel("Normalised G/Gb")
    ylabel("nL [1/min]")
    
    legend
    
end
end