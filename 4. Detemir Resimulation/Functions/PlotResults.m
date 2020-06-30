function [] = PlotResults(P)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P - patient struct

global C

toDateTime = @(mins) P.data.simTime(1) + minutes(mins);

tArray = P.results.tArray;     % Time of results [min]
dtArray = toDateTime(tArray);  % ''              [datetime]

% Set up figure.
patientLabel = sprintf("Patient %d: ", P.patientNum);
F = PanelFigures(3, 3);


%% Glucose
subplot(4, 1, 1)
hold on

[tG, vG] = GetSimTime(P, P.data.G{3});
dtG = P.data.simTime(1) + minutes(tG);
plt = plot(dtG, vG, 'r*');
plt.DisplayName = 'Blood Test';

ppG = griddedInterpolant(tG, vG);
plt = plot(dtArray, ppG(tArray), 'b');
plt.DisplayName = 'Interpolation';

plt = plot(dtArray, P.results.G, 'k');
plt.DisplayName = 'Model Prediction';

lineBounds = ylim;
for ii = 1:length(P.results.nLxLSplits)
    split = toDateTime(P.results.nLxLSplits(ii));
    L = line([split split], lineBounds);
    L.LineWidth = 0.5;
    L.Color = 'k';
    L.HandleVisibility = 'off';
end

title([patientLabel 'Plasma Glucose'])
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')
legend()

datetick('x')
ylim([4 15])


%% Insulin + Detemir
subplot(4, 1, 2)
hold on

[tITotal, vITotal] = GetSimTime(P, P.data.ITotal);  % [pmol/L]
dtITotal = P.data.simTime(1) + minutes(tITotal);
plt = plot(dtITotal, vITotal, 'r*');
plt.DisplayName = 'Blood Test';

ppITotal = griddedInterpolant(tITotal, vITotal);
plt = plot(dtArray, ppITotal(tArray), 'b');
plt.DisplayName = 'Interpolation';

ITotal = C.mU2pmol(P.results.I + P.results.IDF);  % [mU/L] -> [pmol/L]
plt = plot(dtArray, ITotal, 'k');
plt.DisplayName = 'Model Prediction';

lineBounds = ylim;
for ii = 1:length(P.results.nLxLSplits)
    split = toDateTime(P.results.nLxLSplits(ii));
    L = line([split split], lineBounds);
    L.LineWidth = 0.5;
    L.Color = 'k';
    L.HandleVisibility = 'off';
end

title([patientLabel 'Plasma Insulin + Detemir'])
xlabel('Time')
ylabel('Plasma Insulins, I + IDF [pmol/L]')
legend()

datetick('x')


%% Insulin Senstivity
subplot(4, 1, 3)
plot(dtArray, P.results.SI, 'k')

title([patientLabel 'Insulin Sensitivity'])
xlabel('Time')
ylabel('$S_I$ [L/mU/min]')

datetick('x')


%% Insulin Secretion
subplot(4, 1, 4)
plot(dtArray, P.results.Uen, 'k')

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

