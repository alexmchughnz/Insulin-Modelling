set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineLineWidth', 2.0);
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 8 6]);

DATAPATH  = fullfile(pwd, '..', 'Data');
FUNCPATH  = fullfile(pwd, 'Functions');
MODELPATH = fullfile(pwd, 'Models');
RESULTPATH = fullfile(pwd, 'Results');
SCRIPTPATH = fullfile(pwd, 'Scripts');

addpath(genpath(DATAPATH));
addpath(genpath(FUNCPATH));
addpath(genpath(MODELPATH));
addpath(genpath(RESULTPATH));
addpath(genpath(SCRIPTPATH));

PATIENTFILE = 'patient_master.xlsx';

save config DATAPATH FUNCPATH MODELPATH RESULTPATH
disp('Config updated.')
clear