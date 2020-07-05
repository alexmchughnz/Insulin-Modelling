%% main.m (Detemir Resimulation)
% (This code generates a bunch of plots. If you want to hide them, go to
% makedebugplots.m and change any 1s to 0s!)

% Author : Alex McHugh
% Created: 04/06/2020

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
makedata

load config

%% Load Data
patientNums = [1];

% Generate patient data structs.
loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
patients = cell(size(patientNums));
for ii = 1 : length(patientNums)
    n = patientNums(ii);
    patients{ii} = loadpatient(n);
end

%% Calculate Patient/Time Dependent Parameters
for ii = 1:length(patients)
    fprintf("\n\n===PATIENT %d===\n", patients{ii}.patientNum)
    
    % Solve for dependent parameters.    
    patients{ii} = EstimateInsulinSecretion(patients{ii});  % (Uen)
    
%     patients{ii} = FindOptimalHepaticClearance(patients{ii}, 'load', 0);  % (nL, xL) by search
    patients{ii} = FitHepaticClearance(patients{ii}, 'single');  % (nL, xL) by MLR
    
    patients{ii} = FindGutEmptyingRate(patients{ii});       % (d2)
    
    patients{ii} = FitInsulinSensitivity(patients{ii}, true);     % (SI)
    
    % Forward simulate models.
    patients{ii} = SolveSystem(patients{ii});
    
end

%% Plot Results
clear PlotResults PanelFigures

for ii = 1:length(patients)
    P = patients{ii};
    PlotResults(P);
end
PanelDebugPlots();
