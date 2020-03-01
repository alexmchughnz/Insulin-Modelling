%% main.m (Detemir)

% Author : Alex McHugh
% Created: 19/02/2020
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineLineWidth', 2.0);

clear
close all
clc

makeconfig
makedata
makeparameters

load config

%% Load Data
global C GI ID GC SC
load('parameters.mat', 'C', 'GI', 'ID', 'GC', 'SC')   

loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
patients = {loadpatient(1), loadpatient(3), loadpatient(4)};

%% Calculate patient/Time Dependent parameters
for ii = 1:length(patients)
    patients{ii} = EstimateInsulinSecretion(patients{ii});
    patients{ii} = FitInsulinSensitivity(patients{ii});
    patients{ii} = SolveSystem(patients{ii});
end

%% Plot Results
for ii = 1:length(patients)
    P = patients{ii};
    
    %Glucose fit
    range = 2:length(P.G{3}.time) - 1;
    
    figure(1)
    plot(P.G{3}.time(range), P.G{3}.value(range),'r*');
    hold on
    plot(P.results.time, P.results.G, 'k');
    legend('Blood Test','Model')
    title('Plasma Glucose')
    xlabel('Time')
    ylabel('Plasma Glucose, G [mmol/L]')
    datetick('x')
    ylim([4 15])

    %Insulin fit
    ITotal = C.IU2mol*(P.results.I*1e-3)*1e+12 + P.results.IDF;
    figure(2)
    hold on
    plot(P.I.time,P.I.value,'r*')
    plot(P.results.time, ITotal, 'k')
    legend('Blood Test','Model')
    title('Plasma Insulin')
    xlabel('Time')
    ylabel('Plasma Insulin, I [pmol/L]')
    datetick('x')
    % ylim([0 2000])

    %SI fit
    range = 1:2160;
    
    figure(3);
    plot(P.results.time(range), P.SI(range), 'k')
    title('Insulin Sensitivity')
    xlabel('Time')
    ylabel('$S_I$ [L/mU/min]')
    datetick('x')
    % ylim([8e-4 12.5e-4])

    figure(4);
    plot(P.Uen.time, P.Uen.value, 'k')
    title('Estimated Endogenous Insulin Secretion')
    xlabel('Time')
    ylabel('$U_{en}$ [mU/min]')
    datetick('x')
    ylim([0 350])
    
    fprintf("P%d: Plotted results.\n", P.patientNum)
    pause
end

