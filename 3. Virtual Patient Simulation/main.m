%% main.m (Virtual Patient Simulation)
% Monte Carlo Simulation with fast-acting insulin.
% Based on Bekisz, S., Holder-Pearson, L., Chase, J. G., & Desaive, T. (n.d.).,
% In silico validation of a new model-based oral-subcutaneous insulin
% sensitivity testing through monte carlo sensitivity analyses.

% Author : Alex McHugh
% Created: 21/04/2020

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

results = SolveSystem();

%% Plot Results
PlotResults(results);
