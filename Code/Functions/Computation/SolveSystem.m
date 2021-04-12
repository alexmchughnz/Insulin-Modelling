function P = SolveSystem(P, allowPlots)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

DP = DebugPlots().SolveSystem;

if ~exist('allowPlots', 'var')
    allowPlots = false;
end

PrintStatusUpdate(P, "Solving entire system!")

%% Models
P = GIModel(P);

if P.data.IType == "detemir"
    P = IDModel(P);
elseif P.data.IDelivery == "subcutaneous"
    P = SCModel(P);
end

P = GCModel(P);


%% Fit Evaluation
% Insulin
[tI, vI] = GetData(P.data.I);  % [mU/L]
[~, fitI] = GetResultsSample(P, tI, P.results.I);

IMAPE = mean(abs(vI - fitI)./vI);
P.results.fits.insulinMAPE = IMAPE;

ISSE = sum((vI - fitI).^2);
P.results.fits.insulinSSE = ISSE;

% Glucose
[tG, vG] = GetData(P.data.G);  % [mU/L]
[~, fitG] = GetResultsSample(P, tG, P.results.G);

GMAPE = mean(abs(vG - fitG)./vG);
P.results.fits.glucoseMAPE = GMAPE;

GSSE = sum((vG - fitG).^2);
P.results.fits.glucoseSSE = GSSE;


%% Debug Plots
if allowPlots
    tArray = P.results.tArray;
    
    if DP.Glucose
        MakeDebugPlot("Plasma Glucose", P, DP);
        
        [tG, vG] = GetData(P.data.G);
        plt = plot(tG, vG, 'g*');
        plt.DisplayName = 'Plasma Sample';        
       
        plt = plot(tArray, P.results.G, 'k');
        plt.DisplayName = 'Model Prediction';
        
        if isfield(P.data, 'tGBolus')
            plt = line([P.data.tGBolus P.data.tGBolus], ylim, ...
                       'Color', 'g', ...
                       'LineStyle', '--');
            plt.DisplayName = 'Bolus Input';
        end
        
        xlabel('Time')
        ylabel('Plasma Glucose, G [mmol/L]')
        legend()
        
        ylim([4 15])
    end
    
    if DP.Insulin
        MakeDebugPlot("Plasma Insulin", P, DP);
        
        if P.source == "Detemir"
            [tI, vI] = GetData(P.data.ITotal);  % [mU/L]
            I = P.results.I + P.results.IDF;  % [mU/L]
            
            pltylabel = 'Plasma Insulins, I + IDF [mU/L]';            
        else
            [tI, vI] = GetData(P.data.I);  % [mU/L]
            I = P.results.I;  % [mU/L]
            
            pltylabel = 'Plasma Insulin, I [mU/l]';
        end
        
        plt = plot(tI, vI, 'r*');
        plt.DisplayName = 'Plasma Sample';    
        
        plt = plot(tArray, I, 'k');
        plt.DisplayName = 'Model Prediction';   
        
        if isfield(P.data, 'tIBolus')
            x = [P.data.tIBolus'; P.data.tIBolus'];
            y = repmat(ylim', numel(P.data.tIBolus), 1);
            
            plt = line(x, y, ...
                       'Color', 'r', ...
                       'LineStyle', '--');
            plt(1).DisplayName = 'Bolus Input';
            for pp = 2:length(plt)
               plt(pp).HandleVisibility = 'off'; 
            end
        end
        
        xlabel('Time [min]')
        ylabel(pltylabel)
        legend()
    end
end

end
