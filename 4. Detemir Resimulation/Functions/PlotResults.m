function [] = PlotResults(P)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P - patient struct

global C

patientLabel = sprintf("Patient %d: ", P.patientNum);
tArray = P.data.simTime(1) + minutes(P.results.tArray);  % Time of results [datetime]

% Set up figure.
F = PanelFigures(3, 3);


%% Glucose
subplot(4, 1, 1)
hold on

[tG, vG] = GetSimTime(P, P.data.G{3});
tG = P.data.simTime(1) + minutes(tG);
plt = plot(tG, vG, 'r*');
plt.DisplayName = 'Blood Test';

plt = plot(tArray, P.results.G, 'k');
plt.DisplayName = 'Model Prediction';

for ii = 1:length(P.iiSplits)
    split = P.iiSplits(ii);
    L = line([split split], ylim);
    L.LineWidth = 0.5;
    L.Color = 'k';
end

title([patientLabel 'Plasma Glucose'])
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')
legend()

ylim([4 15])


%% Insulin + Detemir
subplot(4, 1, 2)
hold on

[tITotal, vITotal] = GetSimTime(P, P.data.ITotal);  % [pmol/L]
tITotal = P.data.simTime(1) + minutes(tITotal);
plt = plot(tITotal, vITotal, 'r*');
plt.DisplayName = 'Blood Test';


ITotal = C.mU2pmol(P.results.I + P.results.IDF);  % [mU/L] -> [pmol/L]
plt = plot(tArray, ITotal, 'k');
plt.DisplayName = 'Model Prediction';

IInterp = C.mU2pmol(P.ppI(P.results.tArray));
plt = plot(tArray, IInterp, 'b');
plt.DisplayName = 'Interpolation';

for ii = 1:length(P.iiSplits)
    split = P.iiSplits(ii);
    L = line([split split], ylim);
    L.LineWidth = 0.5;
    L.Color = 'k';
end

title([patientLabel 'Plasma Insulin + Detemir'])
xlabel('Time')
ylabel('Plasma Insulins, I + IDF [pmol/L]')
legend()

datetick('x')


%% Insulin Senstivity
subplot(4, 1, 3)
plot(tArray, P.results.SI, 'k')

title([patientLabel 'Insulin Sensitivity'])
xlabel('Time')
ylabel('$S_I$ [L/mU/min]')

datetick('x')


%% Insulin Secretion
subplot(4, 1, 4)
plot(tArray, P.results.Uen, 'k')

title([patientLabel 'Estimated Endogenous Insulin Secretion'])
xlabel('Time')
ylabel('$U_{en}$ [mU/min]')

datetick('x')
ylim([0 300])


%%
path = fullfile("Plots", "patient" + P.patientNum);
savefig(F, path);
fprintf("P%d: Plotted results.\n", P.patientNum)

end

