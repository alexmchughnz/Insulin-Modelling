USERNAME = 'adm181';
DATAPATH = fullfile('C:\Users', USERNAME, 'OneDrive\PhD\Insulin Modelling\1. Detemir\My Code\Data');
FUNCPATH = fullfile('C:\Users', USERNAME, '\OneDrive\PhD\Insulin Modelling\1. Detemir\My Code\Functions');

addpath(DATAPATH);
addpath(FUNCPATH);


PATIENTFILE = 'patient_master.xlsx';

save config DATAPATH FUNCPATH
disp('Config updated.')
clear