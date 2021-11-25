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

Trial.timepoint = tic;

%% Set Up Trial
Trial.label = "";

Trial.recipe = @FixedxLSplineSim;

Trial.source = "CREBRF2021";
Trial.patients = 'all';
% T.patients = 24;


%% Run
Trial = LoadData(Trial);

% Execute on each patient.
runtime = tic;

patientSetIn = Trial.patientSet;
patientSetOut = {};

for ii = 1 : numel(patientSetIn)    
    runtime = PrintTimeRemaining("Main", runtime, ii, Trial.numPatients, Trial.patientSet{ii}, true);
    
    P = patientSetIn{ii};
  
    try
    PArray = Trial.recipe(P);
    catch
        PrintStatusUpdate(P, "Error in simulation. Continuing to next subject.", true);
        continue
    end
    
    SavePatients(Trial, PArray);
    
    patientSetOut = [patientSetOut; PArray(:)];    
end

Trial.resultSet = patientSetOut;

%% Results
PrintTimeTaken(Trial, "Main");

Trial.timepoint = tic;
PrintResults(Trial);
PrintTimeTaken(Trial, "Results");

