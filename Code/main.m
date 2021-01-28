%% main.m
% Author : Alex McHugh
% Created: 26/01/21

clc
clear
close all
fprintf("Running main - press key to start.\n")
% pause
tic

makeconfig
makedebugplots
makeparameters

%% Select Data
% patientNums = [1];
% source = "Detemir";

% patientNums = [24];
% source = "DISST";

% patientNums = [146];
% source = "CREBRF";

patientNums = [5];
source = "OGTTLui";

%% Load Data
patientSet = LoadData(source, patientNums);

%% Run
for ii = 1:length(patientSet)
    patientSet{ii} = SimpleSimulation(patientSet{ii});
end

%% Results
save = true;
PrintResults(patientSet, save);
