%% main.m
% Author : Alex McHugh
% Created: 26/01/21

fprintf("Running main - press key to start.\n")
pause

clc
clear
close all

config

tic
tStart = tic;

%% Select Data
% patientNums = [1];
% source = "Detemir";

% patientNums = [24];
% source = "DISST";

% patientNums = [146];
% source = "CREBRF";

% patientNums = [1 2 4 5 14 16 22 23 25 30];
% patientNums = 25;
patientNums = [5 25 30];
source = "OGTTLui";

%% Load Data
patientSet = LoadData(source, patientNums);
patientSetOut = {};

%% Run
% Select recipe to run on each patient.
recipeFunction = @SimpleSimulation;

% Execute on each patient.
for ii = 1:length(patientSet)
    patientsOut = recipeFunction(patientSet{ii});
    patientSetOut = [patientSetOut patientsOut];
end

%% Results
saveResults = true;
tResults = tic;

PrintTimeTaken("Main", patientSet, tStart);

SavePatients(patientSetOut);
PrintResults(patientSetOut, saveResults);

PrintTimeTaken("Results", patientSet, tResults);

