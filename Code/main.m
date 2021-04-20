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
% recipeFunction = @GridSearchSim;
recipeFunction = @JLKGridSearchSim;
% resultsTag = "Penalised";

% %% Select Data
% patientNums = 'all';
% patientNums = [14];
patientNums = [1 4 14 22 23 25 30];
source = "OGTTLui";

%% Load Data
patientSet = LoadData(source, patientNums, false);
patientSetOut = {};

%% Run
% Execute on each patient.
numPatients = length(patientSet);
runtime = tic;
for ii = 1:numPatients
    runtime = PrintTimeRemaining("Main", runtime, ii, numPatients, patientSet{ii}, true);
    
    % TEMP
    nLArray = [0.1 0.1 0.1...
               0.1 0.1...
               0.15 0.06 0.14 0.06 0.16];
    patientSet{ii}.results.nL = nLArray(ii);
    % \TEMP

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

