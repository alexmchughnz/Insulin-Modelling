
figure();


%% DISST


patientNums = [1 8 5 7 2 3 13 9 10 24];
source = "DISST";
patients = makedata(source, patientNums, false);

PlotProtocol("Trial A", 1, patients)


%% CREBRF
patientNums = [33 79 115 160 169 186 194 196 216 251];  % My chosen 10
source = "CREBRF";
patients = makedata(source, patientNums, false);

PlotProtocol("Trial B", 2, patients)

%% Prettying
subplot(2, 1, 1)
legend()
pause(0.1)

%% Function
function [] = PlotProtocol(plotTable, n, patients)

green = [0.4660 0.6740 0.1880];
blue = [0 0.4470 0.7410];
alpha = 0.3;

subplot(2,1,n)
title(plotTable)
hold on


% Glucose Measurements
for ii = 1:length(patients)
    tG(ii, :) = patients{ii}.data.G.time;
end
center = mean(tG);
err = range(tG);
% plt = errorbar(center, 0.5*ones(size(center)), err, 'horizontal', 'kx');
plt = scatter(center, 0.5*ones(size(center)), 'kx');
plt.DisplayName = "Measurement times";


% Glucose Bolus
tGBolus = cellfun(@(P) P.data.tGBolus, patients);
rectangle('Position', [min(tGBolus) 0 range(tGBolus) 1], ...
    'FaceColor', [green alpha], ...
    'EdgeColor', 'none')
plt = scatter(NaN, NaN, ...
    'MarkerFaceColor', green, ...
    'MarkerFaceAlpha', alpha, ...
    'MarkerEdgeColor', 'none');
plt.DisplayName = "10 g glucose bolus";
% Insulin Bolus
tIBolus = cellfun(@(P) P.data.tIBolus, patients);
rectangle('Position', [min(tIBolus) 0 range(tIBolus) 1], ...
    'FaceColor', [blue alpha], ...
    'EdgeColor', 'none')
plt = scatter(NaN, NaN, ...
    'MarkerFaceColor', blue, ...
    'MarkerFaceAlpha', alpha, ...
    'MarkerEdgeColor', 'none');
plt.DisplayName = "1 U insulin bolus";


xlabel('Time [min]')
xlim([0 40])
set(gca, 'yticklabel', [])
grid
end