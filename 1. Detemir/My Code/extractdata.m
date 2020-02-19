
patientNums = [1 3 4];


% STUB: grab data from previous code structs
for ii = 1:length(patientNums)
    patientNum = patientNums(ii);
    filename = sprintf("sys%d.mat", patientNum);
    load(filename);    
    data = sys.Data;
    
    
    P.patientNum = data.PtNo;
    P.mass = data.pt_mass;
    
    P.CPep.value = data.Cpep;
    P.CPep.time = data.Cpep_time;    
    
    P.G{1}.value = data.bg1;
    P.G{1}.time = data.bg1_time;
    P.G{2}.value = data.bg2;
    P.G{2}.time = data.bg2_time;
    P.G{3}.value = data.bg3;
    P.G{3}.time = data.bg3_time;
    
    P.I.value = data.PlasmaI;
    P.I.time  = data.PlasmaI_time;
    
    P.meal.duration = data.meal_durations;
    P.meal.startTime = data.meal_start;
    P.meal.carbs = data.carbs;
    P.meal.sugar = data.sugar;    
    
    
    filename = sprintf("patient%d.mat", patientNum);
    save(filename, '-struct', 'P');
    clear P
end
%\STUB




% 
% %% Indices
% ROW_1 = 2;
% 
% COL_T = 3;
% COL_G = 4;
% COL_I = 5;
% COL_CPEP = 6;
% COLS_MEAL = [7:14, 16];
% COL_G
% COL_BOLUS = 15;
% 
% 
% 
% 
% %% Load Patient Data
% 
% 
% for ii = 1:length(patientNums)
%     patientNum = patientNums(ii);
%     S = struct('patientNum', patientNum); 
%     [numdata, textdata] = xlsread(PATIENTFILE, patientNum);
%     
%     % Trial Time Data
%     
%     
%     % C-peptide Data
%     
%     % Insulin Bolus Data
%     
%     
%     % Meal Data
%     
%     % Blood Glucose Data
%     
%     S.G = 
%     
%     % Plasma Insulin Data
%     
%     % Fasting Blood Glucose
%     
%     
% end
% 
% 
% %% Load CGM Data