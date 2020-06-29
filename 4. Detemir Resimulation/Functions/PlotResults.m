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
time = P.data.simTime(1) + P.results.tArray/24/60;  % Time of results [datetime]

% Set up figure.
F = PanelFigures(3, 3);

%% Glucose
[tG, vG] = GetSimTime(P, P.data.G{3});

subplot(4, 1, 1)
plot(tG, vG,'r*');
hold on
plot(P.results.tArray, P.results.G, 'k');

for ii = 1:length(P.iiSplits)
    split = P.iiSplits(ii);
    L = line([split split], ylim);
    L.LineWidth = 0.5;
    L.Color = 'k';
end


legend('Blood Test','Model')
title([patientLabel 'Plasma Glucose'])
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')

ylim([4 15])

%% Insulin
ITotal = C.mU2pmol(P.results.I + P.results.IDF);
IInterp = C.mU2pmol(P.ppI(P.results.tArray));
[tI, vI] = GetSimTime(P, P.data.I);

subplot(4, 1, 2)
hold on
plot(tI, vI, 'r*')
plot(P.results.tArray, IInterp, 'b')
plot(P.results.tArray, ITotal, 'k')

for ii = 1:length(P.iiSplits)
    split = P.iiSplits(ii);
    L = line([split split], ylim);
    L.LineWidth = 0.5;
    L.Color = 'k';
end

legend('Blood Test', 'Interpolation', 'Model')

title([patientLabel 'Plasma Insulin'])
xlabel('Time')
ylabel('Plasma Insulin, I [pmol/L]')
% datetick('x')

%% Insulin Senstivity
range = 1 : length(time);

subplot(4, 1, 3)
plot(time(range), P.results.SI(range), 'k')

title([patientLabel 'Insulin Sensitivity'])
xlabel('Time')
ylabel('$S_I$ [L/mU/min]')
datetick('x')

%% Insulin Secretion
subplot(4, 1, 4)
plot(time, P.results.Uen, 'k')

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

