function patientSet = MakeDetemir(patientNums)
% Function for loading Detemir data.
% INPUTS:
%   patientNums - array of patients numbers to load
% OUTPUT:
%   patientSet  - cell array of patient structs 

global CONFIG 
global C SC
global DEBUGPLOTS
if ~exist('allowPlots', 'var')
    allowPlots = true;
end

patientSet = cell(size(patientNums));
source = "Detemir";

for ii = 1:length(patientNums)    
    %% Meta
    patientNum = patientNums(ii);
    filename = sprintf("sys%d.mat", patientNum);
    load(fullfile(CONFIG.DATAPATH, source, filename));
    data = sys.Data;
    
    P.source = source;
    P.patientNum = data.PtNo;
    P.patientCode = sprintf("P%d", P.patientNum);
    P.data.mass = data.pt_mass; % Patient mass [kg]
    
    BMIValues = [27.2 000 27.1 29.1]; %HARDCODED
    P.data.BMI = BMIValues(ii);
    
    ageValues = [73 000 74 75]; %HARDCODED
    P.data.age = ageValues(ii);
    
    P = GetCPeptideParameters(P);
    

    GetMins = @(dt) minutes(dt - sys.sim_start_t);
    
    P.data.trialTime     = [GetMins(sys.trial_start_t), GetMins(sys.trial_end_t)];
    P.data.trialDuration = @() diff(P.data.trialTime);
    
    P.data.simTime     =  [GetMins(sys.sim_start_t), GetMins(sys.sim_end_t)];
    P.data.simDuration =  @() diff(P.data.simTime);
    
    P.results.tArray = (0 : P.data.simDuration())';
    
    
    P.data.CPep.value = data.Cpep;        % C-peptide reading [pmol/L]
    P.data.CPep.time = GetMins(data.Cpep_time);    % Time of C-peptide reading [mins]
    
    P.data.GOther{1}.value = data.bg1;         % Blood glucose reading [mmol/L]
    P.data.GOther{1}.time = data.bg1_time;     % Time of blood glucose reading [mins]
    P.data.GOther{2}.value = data.bg2;
    P.data.GOther{2}.time = data.bg2_time;
    P.data.G.value = data.bg3;
    P.data.G.time = GetMins(data.bg3_time);
    
    P.data.ITotal.value = data.PlasmaI;        % Plasma insulin [pmol/L]
    P.data.ITotal.time  = GetMins(data.PlasmaI_time);
    
    P.data.IBolus = @(~) 0;  % No fast-I bolus here!
    
    vIDBolus = sys.SC.Ibolus;  % Insulin bolus [mU]
    tIDBolus = sys.SC.T;       % Time of bolus delivery [min]
    TIDBolus = 5;              % Period of bolus action [min]
    % Bolus as function of time, value spread over period.
    % Active if time within period.
    P.data.IDBolus = @(t) ((tIDBolus <= t) && (t < tIDBolus+TIDBolus)).*vIDBolus/TIDBolus;
    
    P.data.meal.durations = data.meal_durations;  %[min]
    P.data.meal.startTimes = data.meal_start;     %[datetime]
    P.data.meal.carbs = data.carbs;               %[g]
    P.data.meal.sugar = data.sugar;               %[g]
    
    day1 = sys.sim_start_t;    % Day 1 start, at sim start time [datetime]
    [Y, M, D] = ymd(day1 + 1);
    day2 = datetime(Y, M, D);  % Day 2 start, at midnight [datetime]
    tFast = minutes(day2 - day1); % Time when reading 2 replaces reading 1 [min]
    GFast1 = sys.GC.fasting_bg1;
    GFast2 = sys.GC.fasting_bg2;
    P.data.GFast = @(t) (t < tFast)*GFast1 + (t >= tFast)*GFast2;
    
    P.results.nLxLFitBounds = [];
    
    %% GInfusion Data (for P1)
    P.data.GInfusion = zeros(size(P.results.tArray)); % By default, no infusion.
    
    if (P.patientNum == 1)
        % Information about infusion.
        MAGIC_DEXTROSE_NUMBER = 1.54;  % Assume this is some kind of "how much
        % glucose from 5% dextrose" factor.
        
        duration = 12;          % Duration of infusion [hrs]
        duration = duration*60; % ''                   [min]
        
        startTime  = datetime('31/03/2017 05:15');
        preSimTime = abs(GetMins(startTime)); % How long infusion ran before sim [min]
        
        startTime = P.data.simTime(1) - preSimTime;         % Start of infusion [min]
        endTime   = duration - preSimTime;  % End of infusion [min]
        
        % Return infusion data.
        iiInfusion = (startTime <= P.results.tArray) & (P.results.tArray < endTime); % 1 if infusion active [logical]
        P.data.GInfusion = iiInfusion .* MAGIC_DEXTROSE_NUMBER/C.MGlucose/60;  % Glucose infusion over sim time [mmol/min]
    end    
    
    
    %% Save
    patientSet{ii} = P;
    clear P
end

