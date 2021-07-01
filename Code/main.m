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
% patientNums = [5];
patientNums = 'best';
source = "OGTTLui";

%% Load Data
showLoadPlots = false;
patientSet = LoadData(source, patientNums, showLoadPlots);
patientSetOut = {};

%% Run
% Execute on each patient.
numPatients = length(patientSet);
runtime = tic;
for ii = 1:numPatients
    runtime = PrintTimeRemaining("Main", runtime, ii, numPatients, patientSet{ii}, true);
   
    patientsOut = recipeFunction(patientSet{ii});    
    
    SavePatients(patientsOut);
    
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

