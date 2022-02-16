function CONFIG = config()

set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultLegendLocation', 'southoutside')
set(groot, 'defaultLegendOrientation', 'vertical')
set(groot, 'defaultLineLineWidth', 2.0);
set(groot, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 20 12]);
format shortG
format compact

lastwarn("");
warning('off', 'MATLAB:table:RowsAddedExistingVars');

clear MakeDebugPlot  % Needed to reset persistents in this function.
clear DebugPlots


CONFIG.ENABLELOADERPLOTS = false;
CONFIG.SAVERESULTS = true;
CONFIG.CLOSEALLFIGURES = false;

CONFIG.MAXTIME = 120;
CONFIG.PATIENTSUBNUMBER = @(base, sub) double(string(base) + "00" + string(sub));
CONFIG.PATIENTFILEFORMAT = @(T, P) sprintf("%s_%s", T.source, P.patientCode);
CONFIG.STATUSDEPTH = 2;
CONFIG.HIGHDETAIL = false;

CONFIG.DEFAULTODEOPTIONS = odeset('RelTol', 1e-5, ...
    'MaxStep', 1, ...
    'InitialStep', 0.5);

CONFIG.FUNCPATH  = fullfile(pwd, 'Functions');
CONFIG.LIBPATH = fullfile(pwd, 'Libraries');
CONFIG.MODELPATH = fullfile(pwd, 'Models');
CONFIG.RECIPEPATH = fullfile(pwd, 'Recipes');
CONFIG.TRIALPATH = fullfile(pwd, 'Trials');
CONFIG.DATAPATH  = fullfile(pwd, '..', 'Data');
CONFIG.RESULTPATH = fullfile(pwd, '..', 'Results');
CONFIG.PLOTPATH = fullfile(pwd, '..', 'Plots');
CONFIG.MODULEPATH  = fullfile(pwd, '..', '..', 'Modules');

addpath(genpath(CONFIG.FUNCPATH));
addpath(genpath(CONFIG.LIBPATH));
addpath(genpath(CONFIG.MODELPATH));
addpath(genpath(CONFIG.RECIPEPATH));
addpath(genpath(CONFIG.TRIALPATH));
addpath(genpath(CONFIG.DATAPATH));
addpath(genpath(CONFIG.RESULTPATH));
addpath(genpath(CONFIG.PLOTPATH));
addpath(genpath(CONFIG.MODULEPATH));


disp('Config updated.')

end