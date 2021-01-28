function P = SolveSystem(P, allowPlots)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

global DEBUGPLOTS

if ~exist('allowPlots', 'var')
    allowPlots = false;
end

PrintStatusUpdate(P, "Solving entire system!")

%% Setup
options = odeset('RelTol', 1e-5, ...
    'AbsTol', 1e-4, ...
    'MaxStep', 0.1);


%% Models
P = GIModel(P, options);

if P.data.IType == "detemir"
    P = IDModel(P, options);
elseif P.data.IDelivery == "subcutaneous"
    P = SCModel(P, options);
end

P = GCModel(P, options);


%% Fit Evaluation
[tI, vI] = GetIFromITotal(P); % [mU/L]
iiData = GetTimeIndex(tI, P.results.tArray);
fitI = P.results.I(iiData);

MAPE = mean(abs(vI - fitI)./vI);
P.results.insulinMAPE = MAPE;


%% Debug Plots
if allowPlots
    DP = DEBUGPLOTS.SolveSystem;
    tArray = P.results.tArray;
    
    if DP.Glucose
        MakeDebugPlot("Plasma Glucose", P, DP);
        
        [tG, vG] = GetSimTime(P, P.data.G);
        plt = plot(tG, vG, 'g*');
        plt.DisplayName = 'Plasma Sample';        
       
        plt = plot(tArray, P.results.G, 'k');
        plt.DisplayName = 'Model Prediction';
        
        xlabel('Time')
        ylabel('Plasma Glucose, G [mmol/L]')
        legend()
        
        ylim([4 15])
    end
    
    if DP.Insulin
        MakeDebugPlot("Plasma Insulin", P, DP);
        
        if P.source == "Detemir"
            [tI, vI] = GetSimTime(P, P.data.ITotal);  % [mU/L]
            I = P.results.I + P.results.IDF;  % [mU/L]
            
            pltylabel = 'Plasma Insulins, I + IDF [mU/L]';            
        else
            [tI, vI] = GetSimTime(P, P.data.I);  % [mU/L]
            I = P.results.I;  % [mU/L]
            
            pltylabel = 'Plasma Insulin, I [mU/l]';
        end
        
        plt = plot(tI, vI, 'r*');
        plt.DisplayName = 'Plasma Sample';    
        
        plt = plot(tArray, I, 'k');
        plt.DisplayName = 'Model Prediction';            
        
        xlabel('Time [min]')
        ylabel(pltylabel)
        legend()
    end
end

end
