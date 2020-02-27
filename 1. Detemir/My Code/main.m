%% main.m (Detemir)

% Author : Alex McHugh
% Created: 19/02/2020
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineLineWidth', 2.0);

clear
close all
clc

makeconfig
makedata
makeparameters

load config

%% Load Data
load('parameters.mat')
loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
P = {loadpatient(1), loadpatient(3), loadpatient(4)};


%% Calculate Patient/Time Dependent Parameters
for ii = 1:length(P)
    P{ii} = EstimateInsulinSecretion(GC, P{ii});
    P{ii} = FitInsulinSensitivity(GI, GC, P{ii});
    P{ii} = SolveSystem(GI, ID, GC, SC, P{ii});
end