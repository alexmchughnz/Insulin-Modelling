function P = SolveSystem(P, allowPlots)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

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
if allowPlots
    P = MakePlots(P);
end

end


function P = MakePlots(P)
tag = "SolveSystem";

tArray = P.results.tArray;

%% Glucose Components

P = AddFigure(P, tag, "GlucoseComponents");
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

%% Plasma Glucose
P = AddFigure(P, tag, "PlasmaGlucose");

[tG, vG] = GetData(P.data.G);
if isfield(P.data, "GCV")
    GError = P.data.GCV .* vG;
    plt = errorbar(tG, vG, GError, 'k.', ...
        'MarkerSize', 5, ...
        'MarkerEdgeColor', 'g', ...
        'MarkerFaceColor', 'g');
else
    plt = plot(tG, vG, 'g*');
end


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

%% Plasma Insulin
P = AddFigure(P, tag, "PlasmaInsulin");

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
    for pp = 2:numel(plt)
        plt(pp).HandleVisibility = 'off';
    end
end

xlabel('Time [min]')
ylabel('Plasma Insulin, I [mU/l]')
ylim([0 +Inf])
legend()

%% Interstitial Insulin
P = AddFigure(P, tag, "InterstitialInsulin");

Q = P.results.Q;  % [mU/L]

plt = plot(tArray, Q, 'k');
plt.DisplayName = 'Model Prediction';

xlabel('Time [min]')
ylabel('Interstitial Insulin, Q [mU/l]')
legend()

%% Coefficient Shapes
P = AddFigure(P, tag, "CoefficientShapes");

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

%% nLxL Shapes
P = AddFigure(P, tag, "Clearance");

subplot(2,1,1)
plt = plot(tArray, P.results.xL*ones(size(tArray)), 'r');
plt.DisplayName = "xL";
ylim([0 1])
ylabel("$x_L$")


subplot(2,1,2)

nL = P.results.nL;
if numel(nL) == 1
    nL = nL * ones(size(tArray));
end

plt = plot(tArray, nL, 'b');
plt.DisplayName = "nL";
ylim([0 1])
ylabel("$$n_L$$")

xlabel("Time [min]")

end
