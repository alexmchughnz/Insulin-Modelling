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

patientNums = [1 2 4 5 14 16 22 23 25 30];
source = "OGTTLui";

%% Load Data
patientSet = LoadData(source, patientNums);

%% Run
for ii = 1:length(patientSet)
    forcenLxL = [0.1 0.67];
    patientSet{ii} = SimpleSimulation(patientSet{ii});
end

%% Results
SavePatients(patientSet);

saveResults = true;
PrintResults(patientSet, saveResults);
