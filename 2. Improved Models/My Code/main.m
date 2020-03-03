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
% Global parameter structs - do not mutate!
global C GI ID GC SC NP
load('parameters.mat', 'C', 'GI', 'ID', 'GC', 'SC', 'NP')

loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
patients = {loadpatient(1), loadpatient(3), loadpatient(4)};

%% Calculate Patient/Time Dependent parameters
for ii = 1:length(patients)
    patients{ii} = EstimateInsulinSecretion(patients{ii});
    patients{ii} = FitInsulinSensitivity(patients{ii});
    patients{ii} = SolveSystem(patients{ii});
end

%% Plot Results
plotnum = 1;
P = patients{plotnum};
PlotResults(P);
fprintf("P%d: Plotted results.\n", P.patientNum)

