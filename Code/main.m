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

patientNums = [5];
source = "OGTTLui";

%% Load Data
patientSetIn = LoadData(source, patientNums);
patientSetOut = {};

%% Run
% Select recipe to run on each patient.
recipeFunction = @AdjustUenSimulation;

% Execute on each patient.
for ii = 1:length(patientSetIn)    
    patientsOut = recipeFunction(patientSetIn{ii});
    patientSetOut = [patientSetOut patientsOut];
end

%% Results
SavePatients(patientSetOut);

saveResults = true;
PrintResults(patientSetOut, saveResults);
