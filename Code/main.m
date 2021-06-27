%% main.m
% Author : Alex McHugh
% Created: 26/01/21

fprintf("Running main - press key to start.\n")
% pause

clc
clear
clear functions
close all

config

tStart = tic;

%% Select Recipe
recipeFunction = @SplineSim;
% resultsTag = "Penalised";

%% Select Data
% patientNums = 23;
patientNums = 'best';
% patientNums = [1 4 14 22 23 25 30];
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
    
    SavePatients(patientsOut(:));
    
    patientSetOut = [patientSetOut; patientsOut(:)];    
end

%% Results
clc
tResults = tic;

PrintTimeTaken("Main", patientSet, tStart);

if ~exist("resultsTag", "var")
    resultsTag = "";
end
PrintResults(patientSetOut, recipeFunction, resultsTag);

PrintTimeTaken("Results", patientSet, tResults);

