%% main.m (Detemir Resimulation)
% (This code generates a bunch of plots. If you want to hide them, go to
% makedebugplots.m and change any 1s to 0s!)

% Author : Alex McHugh
% Created: 04/06/2020

tic
clc
close all
fprintf("Running main - press key to start.\n")
% pause 

clear
clear FitHepaticClearance  % to clear persistents in nL/xL plots...
clear global

makedebugplots
makeconfig
makeparameters

load config

SAVERESULTS = false;

%% Load Data
% patientNums = [1 3 4];
% source = "Detemir";

% patientNums = [3 5 7 8 9 13 14 16 24 25];
% source = "DISST";

% patientNums = [12 128 146 160 166 169 171 196 198 216];  % My chosen 10
% source = "CREBRF";

patientNums = [1 2 4 5 14 16 22 23 25 30];
source = "OGTTLui";

patients = makedata(source, patientNums);


%% Calculate Patient/Time Dependent Parameters
for ii = 1:length(patients)
    fprintf("\n\n===PATIENT %s===\n", patients{ii}.patientCode)
    
    %% Determine Uen.
    patients{ii} = EstimateInsulinSecretion(patients{ii});  % (Uen)
%     patients{ii} = AdjustInsulinSecretion(patients{ii}, 'onetooth', 20);
    
    %% Determine nL/xL.
%     patients{ii} = FindOptimalHepaticClearance(patients{ii}, ... 
%         'load');%, 'grid nL[-0.1 0.775]@0.025 xL[0.075 0.95]@0.025');
    
    % Include this parameter to force fit a specific nL xL value.
    forcenLxL = [0.15 0.67];
    patients{ii} = FitHepaticClearance(patients{ii});  % (nL, xL) by MLR

    %% Analyse data variance.
%     stddev = 5/100; 
%     N = 1000;    
%     AnalyseInsulinVariance(patients{ii}, stddev, N);    
    
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
    SaveOpenFigures(source);
    writetable(T, fullfile(RESULTPATH, source+"table.csv"));
end

PanelDebugPlots(1);
disp(T);