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
% Insulin
[tI, vI] = GetIFromITotal(P); % [mU/L]
iiData = GetTimeIndex(tI, P.results.tArray);
fitI = P.results.I(iiData);

IMAPE = mean(abs(vI - fitI)./vI);
P.results.fits.insulinMAPE = IMAPE;

ISSE = sum((vI - fitI).^2);
P.results.fits.insulinSSE = ISSE;

% Glucose
vG = P.data.G.value;
tG = P.data.G.time;
iiData = GetTimeIndex(tG, P.results.tArray);
fitG = P.results.G(iiData);

GMAPE = mean(abs(vG - fitG)./vG);
P.results.fits.glucoseMAPE = GMAPE;

GSSE = sum((vG - fitG).^2);
P.results.fits.glucoseSSE = GSSE;


%% Debug Plots
if allowPlots
    tArray = P.results.tArray;
    
    if DP.Glucose
        MakeDebugPlot("Plasma Glucose", P, DP);
        
        [tG, vG] = GetSimTime(P, P.data.G);
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
