function [] = PlotResults(P)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P - patient struct

global C

ToDateTime = @(mins) P.data.simTime(1) + minutes(mins);

tArray = P.results.tArray;     % Time of results [min]
dtArray = ToDateTime(tArray);  % ''              [datetime]

% Set up figure.
patientLabel = sprintf("Patient %d: ", P.patientNum);
F = PanelFigures(3, 3);


%% Glucose
subplot(4, 1, 1)
hold on

[tG, vG] = GetSimTime(P, P.data.G{3});
dtG = ToDateTime(tG);
plt = plot(dtG, vG, 'r*');
plt.DisplayName = 'Blood Test';

ppG = griddedInterpolant(tG, vG);
plt = plot(dtArray, ppG(tArray), 'b');
plt.LineWidth = 1;
plt.DisplayName = 'Interpolation';

plt = plot(dtArray, P.results.G, 'k');
plt.DisplayName = 'Model Prediction';

lineBounds = ylim;
for ii = 1:length(P.results.nLxLFitBounds)
    split = ToDateTime(P.results.nLxLFitBounds(ii));
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


%% Glucose Error
ax = subplot(4, 1, 2);

iiG = GetTimeIndex(tG, tArray);
simG = P.results.G(iiG);
GError = 100*abs((simG - vG) ./ vG);
plot(dtG, GError, 'r');

title([patientLabel 'Plasma Glucose Error'])
xlabel('Time')
ylabel('Error [\%]')


%% Insulin + Detemir
subplot(4, 1, 3)
hold on

[tITotal, vITotal] = GetSimTime(P, P.data.ITotal);  % [pmol/L]
dtITotal =ToDateTime(tITotal);
plt = plot(dtITotal, vITotal, 'r*');
plt.DisplayName = 'Blood Test';

ppITotal = griddedInterpolant(tITotal, vITotal);
plt = plot(dtArray, ppITotal(tArray), 'b');
plt.LineWidth = 1;
plt.DisplayName = 'Interpolation';

ITotal = C.mU2pmol(P.results.I + P.results.IDF);  % [mU/L] -> [pmol/L]
plt = plot(dtArray, ITotal, 'k');
plt.DisplayName = 'Model Prediction';

lineBounds = ylim;
for ii = 1:length(P.results.nLxLFitBounds)
    split = ToDateTime(P.results.nLxLFitBounds(ii));
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


%% Insulin Error
ax = subplot(4, 1, 4);

iiITotal = GetTimeIndex(tITotal, tArray);
simITotal = ITotal(iiITotal);
ITotalError = 100*abs((simITotal - vITotal) ./ vITotal);
plot(dtITotal, ITotalError, 'r');

title([patientLabel 'Plasma Insulins Error'])
xlabel('Time')
ylabel('Error [\%]')


%%
path = fullfile("Plots", "patient" + P.patientNum);
savefig(F, path);
fprintf("P%d: Plotted results.\n", P.patientNum)

end

