%% main.m (Virtual Patient Simulation)
% Monte Carlo Simulation with fast-acting insulin.
% Based on Bekisz, S., Holder-Pearson, L., Chase, J. G., & Desaive, T. (n.d.).,
% In silico validation of a new model-based oral-subcutaneous insulin
% sensitivity testing through monte carlo sensitivity analyses. 18.

% Author : Alex McHugh
% Created: 30/04/2020

clear
close all
clc

makeconfig
makeparameters

load config

%% Load Data
% Global parameter structs - do not mutate!
global C GI IN GC
load('parameters.mat', 'C', 'GI', 'IN', 'GC')


%% Calculate Patient/Time Dependent parameters
%     patients{ii} = GetGlucoseInfusion(patients{ii});
%     patients{ii} = PDEstimateInsulinSecretion(patients{ii});
%     patients{ii} = FitHepaticClearance(patients{ii});
%     patients{ii} = FitInsulinSensitivity(patients{ii});
%     patients{ii} = SolveSystem(patients{ii});

%% Plot Results
% plotnum = 1;
% P = patients{plotnum};
% PlotResults(P);
% fprintf("P%d: Plotted results.\n", P.patientNum)
