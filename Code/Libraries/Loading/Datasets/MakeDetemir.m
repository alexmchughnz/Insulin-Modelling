function patientSet = MakeDetemir(patientSet)
% Function for loading Detemir data.
% INPUTS:
%   patientSet  - cell array of patient structs
% OUTPUT:
%   patientSet  - updated cell array of patient structs

CONST = Constants();

source = "Detemir";

%% Generate Patients
for ii = 1:length(patientSet)
    P = patientSet{ii};
    
    %% Load Data
    filename = sprintf("sys%d.mat", P.patientNum);
    try
        load(fullfile(CONFIG.DATAPATH, source, filename));
    catch
        assert(false, "Invalid patient number.")  
    end
    data = sys.Data;
    
    %% Patient Info
    ageValues = [73 NaN 74 75];  % HARDCODED
    P.data.age = ageValues(ii);  % [years]
    
    P.data.mass = data.pt_mass;  % [kg]
    
    BMIValues = [27.2 NaN 27.1 29.1];  % HARDCODED
    P.data.BMI = BMIValues(ii);  % [kg/m^2]
    
    %% Trial Times    
    GetMins = @(dt) minutes(dt - sys.sim_start_t);  
    
    P.data.simTime     =  [GetMins(sys.sim_start_t), GetMins(sys.sim_end_t)];
    P.data.simDuration = floor(diff(P.data.simTime));
    P.results.tArray = (P.data.simTime(1) : P.data.simTime(end))';  % [min]
    
    %% Assay Data
    % Glucose Assay
    P.data.G.value = data.bg3;  % [mmol/L]
    P.data.G.time = GetMins(data.bg3_time);  % [min]
    
    day1 = sys.sim_start_t;    % Day 1 start, at sim start time [datetime]
    [Y, M, D] = ymd(day1 + 1);
    day2 = datetime(Y, M, D);  % Day 2 start, at midnight [datetime]
    tFast = minutes(day2 - day1); % Time when reading 2 replaces reading 1 [min]
    GFast1 = sys.GC.fasting_bg1;
    GFast2 = sys.GC.fasting_bg2;
    P.data.GFast = @(t) (t < tFast)*GFast1 + (t >= tFast)*GFast2;
    
    % Insulin Assay    
    P.data.ITotal.value = CONST.pmol2mU(data.PlasmaI);  % [mU/L]
    P.data.ITotal.time  = GetMins(data.PlasmaI_time);  % [min]
    
    % C-peptide Assay
    P.data.CPep.value = data.Cpep;  % [pmol/L]
    P.data.CPep.time = GetMins(data.Cpep_time);  % [min]
    
    %% Trial Inputs
    % Insulin Bolus
    P.data.IType = "detemir";
    P.data.IDelivery = "subcutaneous";
    
    P.data.vIBolus = sys.SC.Ibolus;  % [mU]
    P.data.tIBolus = sys.SC.T;       % [min]
    P.data.TIBolus = 5;              % [min]    
    
    % Glucose Bolus
    P.data.GDelivery = "enteral";
    
    P.data.vGBolus = (data.carbs) / CONST.MGlucose * 1000;  % [mmol]
    P.data.tGBolus = data.meal_start;
    P.data.TGBolus = data.meal_durations;
    
    % Glucose Infusion
    if (P.patientNum == 1)
        MAGIC_DEXTROSE_NUMBER = 1.54;  % Assume "how much glucose from 5% dextrose" factor.
        
        duration = 12;          % Duration of infusion [hrs]
        duration = duration*60; % ''                   [min]
        
        startTime  = datetime('31/03/2017 05:15');
        preSimTime = abs(GetMins(startTime));  % How long infusion ran before sim [min]
        
        startTime = P.data.simTime(1) - preSimTime;  % Start of infusion [min]
        endTime   = duration - preSimTime;  % End of infusion [min]
        
        % Return infusion data.
        iiInfusion = (startTime <= P.results.tArray) & (P.results.tArray < endTime); % 1 if infusion active [logical]
        P.data.GInfusion = iiInfusion .* MAGIC_DEXTROSE_NUMBER/CONST.MGlucose/60;  % [mmol/min]
    else
        P.data.GInfusion = zeros(size(P.results.tArray)); % By default, no infusion.
    end
    
    %% Other
    P = GetCPeptideParameters(P);
    
    %% Save
    patientSet{ii} = P;
    clear P
end
end
