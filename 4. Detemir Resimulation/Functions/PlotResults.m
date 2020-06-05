function [] = PlotResults(P)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P - patient struct

patientLabel = sprintf("Patient %d: ", P.patientNum);
time = P.simTime(1) + P.results.tArray/24/60;  % Time of results [datetime]

% Set up figure.
persistent n;
if (isempty(n))
    n = 0;
end
screensize = get(0,'screensize');
w = screensize(3);
h = screensize(4);
F = figure(P.patientNum);
F.Position = [n/3*w, 0, w/3, 0.9*h];
n = n + 1;

%% Glucose
range = 2 : length(P.G{3}.time) - 1;

subplot(4, 1, 1)
plot(P.G{3}.time(range), P.G{3}.value(range),'r*');
hold on
plot(time, P.results.G, 'k');

legend('Blood Test','Model')
title([patientLabel 'Plasma Glucose'])
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')

datetick('x')
ylim([4 15])

%% Insulin
ITotal = P.results.I + P.results.IDF;

subplot(4, 1, 2)
hold on
plot(P.I.time,P.I.value,'r*')
plot(time, ITotal, 'k')

legend('Blood Test', 'Model')

title([patientLabel 'Plasma Insulin'])
xlabel('Time')
ylabel('Plasma Insulin, I [pmol/L]')
datetick('x')

%% Insulin Senstivity
range = 1 : length(time)-1;

subplot(4, 1, 3)
plot(time(range), P.SI(range), 'k')

title([patientLabel 'Insulin Sensitivity'])
xlabel('Time')
ylabel('$S_I$ [L/mU/min]')
datetick('x')

%% Insulin Secretion
subplot(4, 1, 4)
plot(time(range), P.Uen.value(range), 'k')

title([patientLabel 'Estimated Endogenous Insulin Secretion'])
xlabel('Time')
ylabel('$U_{en}$ [mU/min]')
datetick('x')

%% 
fprintf("P%d: Plotted results.\n", P.patientNum)

end

