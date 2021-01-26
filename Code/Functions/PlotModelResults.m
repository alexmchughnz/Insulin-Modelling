function PlotModelResults(P)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
% INPUTS:
%   P      - patient struct

global DEBUGPLOTS

% Set up figure.
patientLabel = sprintf("%s: ", P.patientCode);
DP = DEBUGPLOTS.PlotResults;

tArray = P.results.tArray;     % Time of results [min]


%% Plasma Glucose
if DP.Glucose
    MakeDebugPlot(P, DP);
    hold on
    
    [tG, vG] = GetSimTime(P, P.data.G);
    plt = plot(tG, vG, 'g*');
    plt.DisplayName = 'Plasma Sample';
    
%     ppG = griddedInterpolant(tG, vG);
%     plt = plot(tArray, ppG(tArray), 'b');
%     plt.LineWidth = 1;
%     plt.DisplayName = 'Interpolation';
%     
%     plt = plot(tArray, P.results.G, 'k');
%     plt.DisplayName = 'Model Prediction';
    
    title(patientLabel + "Plasma Glucose")
    xlabel('Time')
    ylabel('Plasma Glucose, G [mmol/L]')
    legend()
    
    ylim([4 15])
end

%% Insulin (+ Detemir)
if DP.Insulin
    MakeDebugPlot(P, DP);
    hold on
    
    if P.source == "Detemir"
        [tI, vI] = GetSimTime(P, P.data.ITotal);  % [mU/L]
        I = P.results.I + P.results.IDF;  % [mU/L]
        
        plttitle = patientLabel + "Plasma Insulin + Detemir";
        pltylabel = 'Plasma Insulins, I + IDF [mU/L]';
        
        datetick('x')
        
    else
        [tI, vI] = GetSimTime(P, P.data.I);  % [mU/L]
        I = P.results.I;  % [mU/L]
        
        plttitle = patientLabel + "Plasma Insulin";
        pltylabel = 'Plasma Insulin, I [mU/l]';
    end
    
    plt = plot(tArray, I, 'k');
    plt.DisplayName = 'Model Prediction';
    
    plt = plot(tI, vI, 'r*');
    plt.DisplayName = 'Plasma Sample';
    
    % ppI = griddedInterpolant(tI, vI);
    % plt = plot(pltxarray, ppI(tArray), 'b');
    % plt.LineWidth = 1;
    % plt.DisplayName = 'Interpolation';    
    
    title(plttitle)
    xlabel('Time [min]')
    ylabel(pltylabel)
    legend()
end


%%
fprintf("P%d: Plotted results.\n", P.patientNum)

end

