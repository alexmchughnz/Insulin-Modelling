%% main.m (Detemir Resimulation)

% Author : Alex McHugh
% Created: 04/06/2020
clear
close all
clc

makeconfig
makedata
makeparameters

load config

%% Load Data
% Global parameter structs - do not mutate!
global C GC GI ID SC
load('parameters.mat', 'C', 'GC', 'GI', 'ID', 'SC')

% Patient data structs.
loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
patients = {loadpatient(1), loadpatient(3), loadpatient(4)};

%% Calculate Patient/Time Dependent Parameters
for ii = 1:length(patients)
    patients{ii} = GetGlucoseInfusion(patients{ii});
    patients{ii} = EstimateInsulinSecretion(patients{ii});
    patients{ii} = FitInsulinSensitivity(patients{ii});
    patients{ii} = SolveSystem(patients{ii});
end

%% Plot Results
for ii = 1:length(patients)
    P = patients{ii};
    PlotResults(P);
end
