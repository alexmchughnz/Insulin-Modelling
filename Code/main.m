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
Trial.label = "Tony";

Trial.recipe = @SimpleSim;
Trial.source = "OGTTLui";
Trial.patients = 'best';

extraFigures = "SolveSystem.QLocal";

% Trial.patients = [43 45 18  24  26  39  33];


%% Run
Trial = LoadData(Trial);

for ii = 1:numel(extraFigures)
    figname = extraFigures(ii);
    figpath = split(figname, ".");
    Trial.figureList = setfield(Trial.figureList, figpath{:}, true);
end


patientSetIn = Trial.patientSet;
patientSetOut = {};

runtime = tic; 
for ii = 1 : numel(patientSetIn)    
    POut = RunTrial(Trial, ii);  % Can be cell array.

    if ~isempty(POut)
        SavePatients(Trial, POut);
%         MakeFiguresPowerPoint(Trial, POut);
        SaveFigures(Trial, POut);
        patientSetOut = [patientSetOut POut(:)];  % To handle if multiple structs are output.
    end
    
    close all
    
end

Trial.patientSet = patientSetOut;


%% Results
PrintTimeTaken(Trial, "Main");
SaveResults(Trial);
PanelFigures(Trial);

