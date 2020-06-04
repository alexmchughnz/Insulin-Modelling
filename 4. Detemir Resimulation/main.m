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
global C GI ID GC SC
load('parameters.mat', 'C', 'GI', 'ID', 'GC', 'SC')

% Patient data structs.
loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
patients = {loadpatient(1), loadpatient(3), loadpatient(4)};

%% Calculate Patient/Time Dependent parameters
for ii = 1:length(patients)
    patients{ii} = GetGlucoseInfusion(patients{ii});
    patients{ii} = PDEstimateInsulinSecretion(patients{ii});
    patients{ii} = FitHepaticClearance(patients{ii});
    patients{ii} = FitInsulinSensitivity(patients{ii});
    patients{ii} = SolveSystem(patients{ii});
end

%% Plot Results
plotnum = 1;
P = patients{plotnum};
PlotResults(P);
fprintf("P%d: Plotted results.\n", P.patientNum)
