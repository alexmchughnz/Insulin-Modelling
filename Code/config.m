set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineLineWidth', 2.0);
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 8 6]);

global CONFIG
CONFIG.DATAPATH  = fullfile(pwd, '..', 'Data');
CONFIG.FUNCPATH  = fullfile(pwd, 'Functions');
CONFIG.MODELPATH = fullfile(pwd, 'Models');
CONFIG.RESULTPATH = fullfile(pwd, 'Results');
CONFIG.SCRIPTPATH = fullfile(pwd, 'Scripts');

addpath(genpath(CONFIG.DATAPATH));
addpath(genpath(CONFIG.FUNCPATH));
addpath(genpath(CONFIG.MODELPATH));
addpath(genpath(CONFIG.RESULTPATH));
addpath(genpath(CONFIG.SCRIPTPATH));

CONFIG.PATIENTFORMAT = @(source, num) sprintf("P%s%d.mat", source, num);

disp('Config updated.')
clear