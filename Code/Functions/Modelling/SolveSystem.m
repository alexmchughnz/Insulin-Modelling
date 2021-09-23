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
NI = numel(vI);

IMAPE = mean(abs(vI - fitI)./vI);
P.results.fits.insulinMAPE = IMAPE;

ISSE = sum((vI - fitI).^2);
IMSE = ISSE/NI;
P.results.fits.insulinSSE = ISSE;
P.results.fits.insulinMSE = IMSE;

% Glucose
[tG, vG] = GetData(P.data.G);  % [mU/L]
[~, fitG] = GetResultsSample(P, tG, P.results.G);
NG = numel(vI);

GMAPE = mean(abs(vG - fitG)./vG);
P.results.fits.glucoseMAPE = GMAPE;

GSSE = sum((vG - fitG).^2);
GMSE = GSSE/NG;
P.results.fits.glucoseSSE = GSSE;
P.results.fits.glucoseMSE = GMSE;


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


%% Glucose Components
if DP.GlucoseComponents
    GC = P.parameters.GC;
    
    % Copied from GCModelODE.
    G = P.results.G;
    GFast  = P.data.GFast;
    if P.data.GDelivery == "intravenous"
        GInput = P.results.GBolus + P.data.GInfusion;  % [mmol/min]
    else
        GInput = P.data.GInfusion;  % [mmol/min]
    end
    Qtot = P.results.Q;
    Qtot0 = Qtot(1);
    if P.data.IType == "detemir"
        % Include interstitial Detemir.
        Qtot = Qtot + P.results.QDF;   % [mU/L]
        Qtot0 = Qtot0 + P.results.QDF(1);  % [mU/L]
    end
    
    MakeDebugPlot("Glucose Components", P, DP);
    hold on
    
    int = @(x) cumtrapz(tArray, x);
    
    plt = plot(tArray, int(GInput/GC.VG));
    plt.DisplayName = 'Exogenous Glucose';
    plt = plot(tArray, int(P.results.d2/GC.VG * P.results.P2));
    plt.DisplayName = 'Glucose from Gut';
    plt = plot(tArray, int(GC.EGP/GC.VG * ones(size(tArray))));
    plt.DisplayName = 'EGP';
    plt = plot(tArray, int(-GC.pg * (G - GFast)));
    plt.DisplayName = 'NIM Uptake';
    plt = plot(tArray, int(-P.results.SI .* (G.*Qtot - GFast.*Qtot0)./(1 + GC.alphaG*Qtot)));
    plt.DisplayName = 'IM Uptake (with SI)';
    plt = plot(tArray, int(- GC.CNS/GC.VG * ones(size(tArray))));
    plt.DisplayName = 'CNS Uptake';
    
    xlabel('Time [min]')
    ylabel('Glucose Contribution [mmol/L]')
    legend()
    
    lim = max(abs(ylim));
    ylim([-lim lim])
    
    grid on
end

%% Glucose
if DP.Glucose
    MakeDebugPlot("Plasma Glucose", P, DP);
    
    [tG, vG] = GetData(P.data.G);
    GError = P.data.GCV .* vG;
    plt = errorbar(tG, vG, GError, 'k.', ...
        'MarkerSize', 5, ...
        'MarkerEdgeColor', 'g', ...
        'MarkerFaceColor', 'g');
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

%% Plasma Insulin
if DP.PlasmaInsulin
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

%% Interstitial Insulin
if DP.InterstitialInsulin
    MakeDebugPlot("Interstitial Insulin", P, DP);
    
    Q = P.results.Q;  % [mU/L]
        
    plt = plot(tArray, Q, 'k');
    plt.DisplayName = 'Model Prediction';
    
    xlabel('Time [min]')
    ylabel('Interstitial Insulin, Q [mU/l]')
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
    
    plt = plot(tArray, P.results.Uex(P)/GC.VI, 'b');
    plt.DisplayName = 'Exogenous';
    
    xlabel("Time [min]")
    ylabel("Insulin [mU/L]")
    legend()
end

end
