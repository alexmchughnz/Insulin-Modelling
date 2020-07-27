%% main.m (Detemir Resimulation)
% (This code generates a bunch of plots. If you want to hide them, go to
% makedebugplots.m and change any 1s to 0s!)

% Author : Alex McHugh
% Created: 04/06/2020

tic
clc
close all
fprintf("Running main - press key to start.\n")
pause 

clear
clear FitHepaticClearance  % to clear persistents in nL/xL plots...
clear global

makedebugplots
makeconfig
makeparameters

load config

%% Load Data
patientNums = [1];
source = "DISST";
patients = makedata(source, patientNums);


%% Calculate Patient/Time Dependent Parameters
for ii = 1:length(patients)
    fprintf("\n\n===PATIENT %s===\n", patients{ii}.patientCode)
    
    %% Determine Uen.
    patients{ii} = EstimateInsulinSecretion(patients{ii});  % (Uen)
    
    %% Determine nL/xL.
%     patients{ii} = FindOptimalHepaticClearance(patients{ii}, ... 
%         'grid', [0 0.3], [0.6 1], 0.02);  % (nL, xL) by search
    
    patients{ii} = FitHepaticClearance(patients{ii}, 'single');  % (nL, xL) by MLR

    %% Analyse data variance.
%     stddev = 5/100; 
%     N = 1000;    
%     AnalyseInsulinVariance(patients{ii}, stddev, N);    
    
    %% Find other dependent parameters. 
    patients{ii} = FindGutEmptyingRate(patients{ii});  % (d2)
    
    patients{ii} = FitInsulinSensitivity(patients{ii}, true);  % (SI)
    
    %% Forward simulate models.
    patients{ii} = SolveSystem(patients{ii});
    
end

%% Plot Results
clear PlotResults PanelFigures

for ii = 1:length(patients)
    P = patients{ii};
    PlotResults(P);
end
PanelDebugPlots(1);
