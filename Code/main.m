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

%% Select Recipe
recipeFunction = @SimpleSim;


%% Select Data
% patientNums = 'all';
patientNums = 23;
source = "MockOGTTLui";

%% Load Data
patientSet = LoadData(source, patientNums, false);
patientSetOut = {};

%% Run
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
PrintResults(patientSetOut, recipeFunction);

PrintTimeTaken("Results", patientSet, tResults);

