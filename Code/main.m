%% main.m
% Author : Alex McHugh
% Created: 26/01/21

fprintf("Running main - press key to start.\n")
pause

clc
clear
close all

config

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
patientNums = [25];
source = "OGTTLui";

%% Load Data
patientSet = LoadData(source, patientNums);
patientSetOut = {};

%% Run
% Select recipe to run on each patient.
recipeFunction = @MatchIInputSim;

% Execute on each patient.
numPatients = length(patientSet);
runtime = tic;
for ii = 1:numPatients
    runtime = PrintTimeRemaining("Main", runtime, ii, numPatients, patientSet{ii}, true);

    patientsOut = recipeFunction(patientSet{ii});
    patientSetOut = [patientSetOut patientsOut];
    
end

%% Results
clc
tResults = tic;

PrintTimeTaken("Main", patientSet, tStart);

SavePatients(patientSetOut);
PrintResults(patientSetOut);

PrintTimeTaken("Results", patientSet, tResults);

