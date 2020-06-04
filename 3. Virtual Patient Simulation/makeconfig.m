set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultLineLineWidth', 2.0);

FUNCPATH  = fullfile(pwd, 'Functions');
MODELPATH = fullfile(pwd, 'Models');

addpath(FUNCPATH);
addpath(MODELPATH);

save config FUNCPATH MODELPATH
disp('Config updated.')
clear