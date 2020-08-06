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

SAVERESULTS = false;

%% Load Data
% patientNums = [1 3 4];
% source = "Detemir";

patientNums = [1 8 5 7 2 3 13 9 10 24];
% patientNums = [1];
source = "DISST";

patients = makedata(source, patientNums);


%% Calculate Patient/Time Dependent Parameters
for ii = 1:length(patients)
    fprintf("\n\n===PATIENT %s===\n", patients{ii}.patientCode)
    
    %% Determine Uen.
    patients{ii} = EstimateInsulinSecretion(patients{ii});  % (Uen)
    
    %% Determine nL/xL.
    patients{ii} = FindOptimalHepaticClearance(patients{ii}, ... 
        'load');  % (nL, xL) by search
    
    nL = patients{ii}.results.nL(1);
    xL = patients{ii}.results.xL(1);    
    patients{ii} = FitHepaticClearance(patients{ii}, [nL xL]);  % (nL, xL) by MLR

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

[~, I] = sort(patientNums);
for ii = 1:length(patientNums)
   idx = find(I == ii);
   P = patients{idx};
   fprintf("%s\t %.2f/%.2f\n", P.patientCode, P.results.nL(1), P.results.xL(1));
end

%% Plot Results
clear PlotResults PanelFigures

T = table;

for ii = 1:length(patients)
    P = patients{ii};
    
    PlotResults(P);
    T = TabulateResults(T, P);
end

PanelDebugPlots(1);
disp(T);

if (SAVERESULTS)
    saveopenfigures;
    writetable(T, fullfile(RESULTPATH, 'table.txt'));
end