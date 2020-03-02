% Adapted from "PtDataRead.m".

load config

patientNums = [1 3 4];


% STUB: grab data from previous code structs
for ii = 1:length(patientNums)
    patientNum = patientNums(ii);
    filename = sprintf("sys%d.mat", patientNum);
    load(filename);    
    data = sys.Data;
    
    clear P
    P.patientNum = data.PtNo;
    P.mass = data.pt_mass;           % Patient mass [kg]
    
    P.trialTime     = [sys.trial_start_t, sys.trial_end_t];
    P.trialDuration = sys.trial_max_t;
    
    %NOTE: Unsure what "sim time" is, but it's used here...?
    P.simTime     =  [sys.sim_start_t, sys.sim_end_t];
    P.simDuration =  sys.sim_max_t;
    
    P.CPep.value = data.Cpep;        % C-peptide reading [pmol/L]
    P.CPep.time = data.Cpep_time;    % Time of C-peptide reading [datetime]
    
    P.G{1}.value = data.bg1;         % Blood glucose reading [mmol/L?]
    P.G{1}.time = data.bg1_time;     % Time of blood glucose reading [datetime]
    P.G{2}.value = data.bg2;
    P.G{2}.time = data.bg2_time;
    P.G{3}.value = data.bg3;
    P.G{3}.time = data.bg3_time;
    
    P.I.value = data.PlasmaI;        % Plasma insulin [?]
    P.I.time  = data.PlasmaI_time;
    
    P.IBolus.value = sys.SC.Ibolus;  % Insulin bolus [mU]
    P.IBolus.time  = sys.SC.T; 
    
    P.meal.durations = data.meal_durations;  %[min]
    P.meal.startTimes = data.meal_start;     %[datetime]
    P.meal.carbs = data.carbs;               %[g]
    P.meal.sugar = data.sugar;               %[g]
    
    P.GFast{1} = sys.GC.fasting_bg1;
    P.GFast{2} = sys.GC.fasting_bg2;    
    
    filename = sprintf("patient%d.mat", patientNum);
    save(fullfile(DATAPATH, filename), '-struct', 'P');
    fprintf('P%d: Saved patient data.\n', patientNum);
end
%\STUB

clear




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