%% main.m
% (This code generates a bunch of plots. If you want to hide them, go to
% makedebugplots.m and change any 1s to 0s!)

% Author : Alex McHugh
% Created: 26/01/21

clc
clear
close all
fprintf("Running main - press key to start.\n")
% pause
tic

debugplots
config
parameters

%% Load Data
patientNums = [1 3 4];
source = "Detemir";


patientSet = LoadData(source, patientNums);




pause