%% main.m
% Author : Alex McHugh
% Created: 26/01/21

fprintf("Running main - press key to start.\n")
% pause

clc
clear
close all

Trial.Config = config();
Trial.startTime = tic;

%% Set Up Trial
% Trial.label = "";

Trial.recipe = @SplineSim;
Trial.source = "CREBRF2021";
% Trial.patients = "all";
% Trial.patients = [43 45 18  24  26  39  33];
Trial.patients = 45;


%% Run
Trial = LoadData(Trial);

patientSetIn = Trial.patientSet;
patientSetOut = {};

runtime = tic; 
for ii = 1 : numel(patientSetIn)    
    POut = RunTrial(Trial, ii);

    if ~isempty(POut)
        SavePatients(Trial, POut);
        MakeFiguresPowerPoint(Trial, POut);
        SaveFigures(Trial, POut);
        patientSetOut = [patientSetOut POut(:)];  % To handle if multiple structs are output.
    end
end

Trial.patientSet = patientSetOut;


%% Results
PrintTimeTaken(Trial, "Main");
SaveResults(Trial);
PanelFigures(Trial);

