%% main.m
% Author : Alex McHugh
% Created: 26/01/21

fprintf("Running main - press key to start.\n")
% pause

clc
clear
close all

config

tStart = tic;

%% Select Recipe
% recipeFunction = @SimpleSim;
recipeFunction = @GridSearchSim;
% recipeFunction = @IterateParametersSim;
% resultsTag = "Penalised";

% %% Select Data
patientNums = 'best';
% patientNums = [2 4 5 14 16 22 23 25 30];
% patientNums = [1 4 14 22 23 25 30];
source = "DISST";

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
if ~exist("resultsTag", "var")
    resultsTag = "";
end
PrintResults(patientSetOut, recipeFunction, resultsTag);

PrintTimeTaken("Results", patientSet, tResults);

