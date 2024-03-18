%% main.m
% Author : Alex McHugh

clc
clear
close all

config

%% Select Data
patientNums = [1];
source = "template";

%% Load Data
patientSet = LoadData(source, patientNums);
patientSetOut = {};

%% Run
% Execute on each patient.
numPatients = length(patientSet);

for ii = 1:numPatients
    patientsOut = SimpleSim(patientSet{ii});    
    patientSetOut = [patientSetOut; patientsOut(:)];    
end

%% Results

PrintResults(patientSetOut, @SimpleSim);
