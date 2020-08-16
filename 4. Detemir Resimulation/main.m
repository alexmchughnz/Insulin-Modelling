%% main.m (Detemir Resimulation)
% (This code generates a bunch of plots. If you want to hide them, go to
% makedebugplots.m and change any 1s to 0s!)

% Author : Alex McHugh
% Created: 04/06/2020

tic
clc
close all
fprintf("Running main - press key to start.\n")
pause 

clear
clear FitHepaticClearance  % to clear persistents in nL/xL plots...
clear global

makedebugplots
makeconfig
makeparameters

load config

SAVERESULTS = true;

%% Load Data
% patientNums = [1 3 4];
% source = "Detemir";

% patientNums = [1 8 5 7 2 3 13 9 10 24];
% source = "DISST";

% patientNums = [33 79 115 153 169 186 194 196 216 251];  %Jen O's chosen 10
patientNums = [33 79 147 160 169 186 194 196 216 251];  % My chosen 10
source = "CREBRF";

patients = makedata(source, patientNums);


%% Calculate Patient/Time Dependent Parameters
for ii = 1:length(patients)
    fprintf("\n\n===PATIENT %s===\n", patients{ii}.patientCode)
    
    %% Determine Uen.
    patients{ii} = EstimateInsulinSecretion(patients{ii});  % (Uen)
    
    %% Determine nL/xL.
%     patients{ii} = FindOptimalHepaticClearance(patients{ii}, ... 
%         'load');  % (nL, xL) by search
    
% %     Include this parameter to force best grid search result.
%     forcenLxL = [0.05 0.67];
%     forcenLxL = [patients{ii}.results.nL(1) patients{ii}.results.xL(1)];
    patients{ii} = FitHepaticClearance(patients{ii});  % (nL, xL) by MLR

    %% Analyse data variance.
%     stddev = 5/100; 
%     N = 1000;    
%     AnalyseInsulinVariance(patients{ii}, stddev, N);    
%     
    %% Find other dependent parameters. 
    patients{ii} = FindGutEmptyingRate(patients{ii});  % (d2)
    
    patients{ii} = FitInsulinSensitivity(patients{ii}, true);  % (SI)
    
    %% Forward simulate models.
    patients{ii} = SolveSystem(patients{ii});
    
end


%% Plot Results
clear PlotResults PanelFigures

T = table;

for ii = 1:length(patients)
    P = patients{ii};
    
    PlotResults(P);
    T = TabulateResults(T, P);
end

if (SAVERESULTS)
    saveopenfigures;
    writetable(T, fullfile(RESULTPATH, 'table.txt'));
end

PanelDebugPlots(2);
disp(T);