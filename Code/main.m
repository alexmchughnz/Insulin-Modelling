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
% recipeFunction = @MatchIInputSim;
recipeFunction = @SimpleSim;
% resultsTag = "Normalised";

%% Select Data
% patientNums = 'all';
patientNums = 23;
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

