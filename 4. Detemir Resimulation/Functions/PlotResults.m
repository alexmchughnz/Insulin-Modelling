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
time = P.simTime(1) + P.results.tArray/24/60;  % Time of results [datetime]

% Set up figure.
F = PanelFigures(3, 3);

%% Glucose
range = 2 : length(P.data.G{3}.time) - 1;

subplot(4, 1, 1)
plot(P.data.G{3}.time(range), P.data.G{3}.value(range),'r*');
hold on
plot(time, P.results.G, 'k');

legend('Blood Test','Model')
title([patientLabel 'Plasma Glucose'])
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')

datetick('x')
ylim([4 15])

%% Insulin
ITotal = C.mU2pmol(P.results.I + P.results.IDF);

subplot(4, 1, 2)
hold on
plot(P.data.I.time,P.data.I.value,'r*')
plot(time, ITotal, 'k')

legend('Blood Test', 'Model')

title([patientLabel 'Plasma Insulin'])
xlabel('Time')
ylabel('Plasma Insulin, I [pmol/L]')
datetick('x')

%% Insulin Senstivity
range = 1 : length(time);

subplot(4, 1, 3)
plot(time(range), P.SI(range), 'k')

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
fprintf(" d2 = %.4f", P.d2)
fprintf(" nL = %.4f | %.4f \n",  P.nL(1), P.nL(end))
fprintf(" xL = %.4f | %.4f \n",  P.xL(1), P.xL(end))

end

