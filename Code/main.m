%% main.m
% Author : Alex McHugh
% Created: 26/01/21

fprintf("Running main - press key to start.\n")
% pause

clc
clear
clear functions
close all

Trial.Config = config();
Trial.timepoint = tic;

%% Set Up Trial
Trial.label = "";

Trial.recipe = @SplineSim;

Trial.source = "CREBRF2021";
Trial.patients = 'all';
% Trial.patients = 13;


%% Run
Trial = LoadData(Trial);

% Execute on each patient.
runtime = tic;

patientSetIn = Trial.patientSet;
patientSetOut = {};

for ii = 1 : numel(patientSetIn)    
    runtime = PrintTimeRemaining("Main", runtime, ii, Trial.numPatients, Trial.patientSet{ii}, true);
    
    PIn = patientSetIn{ii};
    
    POut = RunTrial(Trial, PIn);
    if ~isempty(POut)
        SavePatients(Trial, POut);
        patientSetOut = [patientSetOut; POut(:)];  
    end
end

Trial.resultSet = patientSetOut;


%% Results
PrintTimeTaken(Trial, "Main");

Trial.timepoint = tic;
PrintResults(Trial);
PrintTimeTaken(Trial, "Results");

