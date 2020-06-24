%% main.m (Detemir Resimulation)

% Author : Alex McHugh
% Created: 04/06/2020
clear
clear global
close all
clc

global DEBUGPLOT
DEBUGPLOT = false;

makeconfig
makeparameters
makedata

load config

%% Load Data
% Patient data structs.
loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
patients = {loadpatient(1), loadpatient(3), loadpatient(4)};

%% Calculate Patient/Time Dependent Parameters
for ii = 1:length(patients)
    
    % Solve for dependent parameters.    
    patients{ii} = EstimateInsulinSecretion(patients{ii});  
    patients{ii} = FitHepaticClearance(patients{ii});
%     patients{ii}.nL = 0.15;
%     patients{ii}.xL = 0.67;
    
    patients{ii} = GridSearch(patients{ii});
    
    patients{ii} = FitInsulinSensitivity(patients{ii});
    
    % Forward simulate models.
    patients{ii} = SolveSystem(patients{ii});
end

%% Plot Results
clear PlotResults PanelFigures
for ii = 1:length(patients)
    P = patients{ii};
    PlotResults(P);
end
