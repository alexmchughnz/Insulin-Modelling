function [num, txt] = get_patient_data(sheet_num)
%Returns num and txt data for a specific patient from the main patient
%database. Directory at \\file\UsersJ$\jkl59\Home\My Documents\408
%FYP\patient data 
%This will need to be changed for your individual directory.

global USERNAME
% (ADM 18/02/20)
% Changed for my directory.
filename = 'patient_master.xlsx';
fullpath = fullfile('C:\Users', USERNAME, 'OneDrive\PhD\Insulin Modelling\1. Detemir\Original Code');
addpath(fullpath);

[num, txt, ~] = xlsread(filename, sheet_num);

end