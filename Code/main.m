%% main.m
% Author : Alex McHugh
% Created: 26/01/21

clc
clear
close all
fprintf("Running main - press key to start.\n")
% pause
tic

config

%% Select Data
% patientNums = [1];
% source = "Detemir";

% patientNums = [24];
% source = "DISST";

% patientNums = [146];
% source = "CREBRF";

patientNums = [1];
source = "OGTTLui";

%% Load Data
patientSet = LoadData(source, patientNums);

%% Run
for ii = 1:length(patientSet)
    forcenLxL = [0.1 0.67];
    patientSet{ii} = SimpleSimulationForced(patientSet{ii}, forcenLxL);
end

%% Results
SavePatients(patientSet);

saveResults = true;
PrintResults(patientSet, saveResults);
