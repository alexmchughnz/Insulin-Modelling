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

T.timepoint = tic;

%% Set Up Trial
T.label = "";

T.recipe = @FixedxLSplineSim;

T.source = "CREBRF2021";
T.patients = 'all';
% T.patients = 24;

%% Load Data
T = LoadData(T);

%% Run
% Execute on each patient.
runtime = tic;

patientSetIn = T.patientSet;
patientSetOut = {};

for ii = 1 : numel(patientSetIn)    
    runtime = PrintTimeRemaining("Main", runtime, ii, T.numPatients, T.patientSet{ii}, true);
    
    P = patientSetIn{ii};
  
    try
    PArray = T.recipe(P);
    catch
        PrintStatusUpdate(P, "Error in simulation. Continuing to next subject.", true);
        continue
    end
    
    SavePatients(T, PArray);
    
    patientSetOut = [patientSetOut; PArray(:)];    
end

T.resultSet = patientSetOut;

%% Results
PrintTimeTaken(T, "Main");

T.timepoint = tic;
PrintResults(T);
PrintTimeTaken(T, "Results");

