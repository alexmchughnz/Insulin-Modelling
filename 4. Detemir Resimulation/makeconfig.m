set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineLineWidth', 2.0);

DATAPATH  = fullfile(pwd, 'Data');
FUNCPATH  = fullfile(pwd, 'Functions');
MODELPATH = fullfile(pwd, 'Models');

addpath(DATAPATH);
addpath(FUNCPATH);
addpath(MODELPATH);

PATIENTFILE = 'patient_master.xlsx';

save config DATAPATH FUNCPATH MODELPATH
disp('Config updated.')
clear