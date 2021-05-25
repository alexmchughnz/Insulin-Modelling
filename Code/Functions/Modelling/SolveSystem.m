function P = SolveSystem(P, allowPlots)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

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


%% Plotting
if ~exist('allowPlots', 'var')
    allowPlots = false;
end

if allowPlots
    MakePlots(P);
end

end


function MakePlots(P)
DP = DebugPlots().SolveSystem;

tArray = P.results.tArray;

%% Glucose
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

%% Insulin
if DP.Insulin
    MakeDebugPlot("Plasma Insulin", P, DP);
    
    [tI, vI] = GetData(P.data.I);  % [mU/L]
    I = P.results.I;  % [mU/L]
    
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
    ylabel('Plasma Insulin, I [mU/l]')
    legend()
end

%% Coefficient Shapes
if DP.CoefficientShapes
    GC = P.parameters.GC;
    
    MakeDebugPlot("Coefficient Shapes", P, DP);    
    
    [tI, vI] = GetData(P.data.I);  % [mU/L]
    plt = plot(tI, vI, 'r*');
    plt.DisplayName = 'Plasma Sample';
    
    plt = plot(tArray, P.results.Uen/GC.VI, 'g');
    plt.DisplayName = 'Endogenous';
    
    plt = plot(tArray, P.results.IInput/GC.VI, 'b');
    plt.DisplayName = 'Exogenous';
    
    xlabel("Time [min]")
    ylabel("Insulin [mU/L]")
    legend()
end

end
