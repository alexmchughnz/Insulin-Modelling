function [num, txt] = get_patient_data(sheet_num)
%Returns num and txt data for a specific patient from the main patient
%database. Directory at \\file\UsersJ$\jkl59\Home\My Documents\408
%FYP\patient data 
%This will need to be changed for your individual directory.


% (ADM 18/02/20)
% Changed for my directory.
filename = 'patient_master.xlsx';
addpath 'C:\Users\adm181\OneDrive\PhD\Insulin Modelling\1. Detemir\Original Code'

[num, txt, ~] = xlsread(filename, sheet_num);

end