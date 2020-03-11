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


%% Glucose
range = 2 : length(P.G{3}.time) - 1;

figure(1)
plot(P.G{3}.time(range), P.G{3}.value(range),'r*');
hold on
plot(P.results.time, P.results.G, 'k');

legend('Blood Test','Model')
title([patientLabel 'Plasma Glucose'])
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')

datetick('x')
ylim([4 15])

%% Insulin
ITotal = P.results.I + P.results.IDF;

figure(2)
hold on
plot(P.I.time,P.I.value,'r*')
plot(P.results.time, ITotal, 'k')

legend('Blood Test','Model')

title([patientLabel 'Plasma Insulin'])
xlabel('Time')
ylabel('Plasma Insulin, I [pmol/L]')
datetick('x')

%% Insulin Senstivity
range = 1 : length(P.results.time)-1;

figure(3);
plot(P.results.time(range), P.SI(range), 'k')

title([patientLabel 'Insulin Sensitivity'])
xlabel('Time')
ylabel('$S_I$ [L/mU/min]')
datetick('x')

%% Insulin Secretion
figure(4);
plot(P.results.time(range), P.Uen.value(range), 'k')

title([patientLabel 'Estimated Endogenous Insulin Secretion'])
xlabel('Time')
ylabel('$U_{en}$ [mU/min]')
datetick('x')
ylim([0 350])

end

