set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultLegendLocation', 'southoutside')
set(groot, 'defaultLegendOrientation', 'vertical')
set(groot, 'defaultLineLineWidth', 2.0);
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 20 12]);
format shortG
format compact
warning('off', 'MATLAB:table:RowsAddedExistingVars');

clear MakeDebugPlot  % Needed to reset persistents in this function.
clear DebugPlots

global CONFIG
CONFIG.FUNCPATH  = fullfile(pwd, 'Functions');
CONFIG.LIBPATH = fullfile(pwd, 'Libraries');
CONFIG.MODELPATH = fullfile(pwd, 'Models');
CONFIG.RECIPEPATH = fullfile(pwd, 'Recipes');
CONFIG.DATAPATH  = fullfile(pwd, '..', 'Data');
CONFIG.RESULTPATH = fullfile(pwd, '..', 'Results');
CONFIG.PLOTPATH = fullfile(pwd, '..', 'Plots');

addpath(genpath(CONFIG.FUNCPATH));
addpath(genpath(CONFIG.LIBPATH));
addpath(genpath(CONFIG.MODELPATH));
addpath(genpath(CONFIG.RECIPEPATH));
addpath(genpath(CONFIG.DATAPATH));
addpath(genpath(CONFIG.RESULTPATH));
addpath(genpath(CONFIG.PLOTPATH));

CONFIG.PATIENTFORMAT = @(P) sprintf("P%s%d", P.source, P.patientNum);
CONFIG.STATUSDEPTH = 3;
CONFIG.SAVERESULTS = true;
CONFIG.HIGHDETAIL = false;

defaultOptions = odeset('RelTol', 1e-5, ...
    'AbsTol', 1e-4, ...
    'MaxStep', 0.5, ...
    'InitialStep', 0.1);
CONFIG.DEFAULTODEOPTIONS = defaultOptions;

disp('Config updated.')
clear