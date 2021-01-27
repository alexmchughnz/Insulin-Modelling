function patientSet = MakeDISST(patientSet)
% Function for loading DISST data.
% INPUTS:
%   patientSet  - cell array of patient structs
% OUTPUT:
%   patientSet  - updated cell array of patient structs


global CONFIG
global C

source = "DISST";

%% Load Data
opts = spreadsheetImportOptions(...
    'NumVariables', 9, ...
    'DataRange', 'C3:K53', ...
    'VariableNamesRange', 'C2:K2', ...
    'RowNamesRange', 'B3:B53');
opts = setvartype(opts, 'double');
infoTable = readtable(fullfile(CONFIG.DATAPATH, source, "database recent.xls"), opts, ...
    'ReadRowNames', true, ...
    'ReadVariableNames', true);

opts = spreadsheetImportOptions(...
    'NumVariables', 25, ...
    'DataRange', 'A3:Y52', ...
    'VariableNamesRange', '2:2');
opts = setvartype(opts, 'double');
dataTable = readtable(fullfile(CONFIG.DATAPATH, source, "DIST recent.xls"), opts, ...
    'ReadRowNames', true, ...
    'ReadVariableNames', true);

nMeas = 5;  % Number of measurements.
    
%% Generate Patients
for ii = 1:length(patientSet)
    P = patientSet{ii};
    code = dataTable.Properties.RowNames{ii};
    assert(P.patientNum <= 50, "Invalid patient number.")
    
    %% Patient Info
    P.data.age = infoTable{code, "age_years_"};  % [years]
    P.data.mass = infoTable{code, "weight_kg_"};  % [kg]
    P.data.height = infoTable{code, "height_cm_"};  % [cm]
    P.data.BMI = infoTable{code, "bmi"};  % [kg/m^2]
    
    %% Trial Times
    measTimes = dataTable{code, repmat("time", 1, nMeas) + (1:nMeas)}'/60;  % Time of measurement [min]
    P.data.simTime = [floor(min(measTimes)), ceil(max(measTimes))];
    P.data.simDuration =  @() floor(diff(P.data.simTime));
    P.results.tArray = (P.data.simTime(1) : 1/60 : P.data.simTime(end))';
    
    %% Assay Data
    % Glucose Assay
    P.data.G.value = dataTable{code, repmat("G", 1, nMeas) + (1:nMeas)}';  % [mmol/L]
    P.data.G.time = measTimes;  % [min]
    P.data.GFast = @(~) P.data.G.value(1);  % [mmol/L]
    
    % Insulin Assay
    P.data.I.value = dataTable{code, repmat("I", 1, nMeas) + (1:nMeas)}';  % [mU/L]
    P.data.I.time = measTimes;  % [min]
    
    % C-peptide Assay
    P.data.CPep.value = dataTable{code, repmat("C", 1, nMeas) + (1:nMeas)}';  % [pmol/L]
    P.data.CPep.time = measTimes;  % [min]
    
    %% Trial Inputs
    MakeBolusFunction = @(value, time, period) (@(tArray) arrayfun(@(t) dot(value/period, (time<=t & t<(time+period))), tArray));
    
    % Insulin Bolus
    vIBolus = dataTable{code, "IB"} * 1e+3;  % [mU]
    tIBolus = dataTable{code, "timeIB"}/60;  %  [min]
    TIBolus = 1;  % [min]
    P.data.IBolus = MakeBolusFunction(vIBolus, tIBolus, TIBolus);  % [mU/min]

    % Glucose Bolus
    vGBolus = dataTable{code, "GB"} / C.MGlucose * 1e+3;  % [mmol]
    tGBolus = dataTable{code, "timeGB"}/60;  % [min]
    TGBolus = 1;  % [min]
    P.data.GBolus = MakeBolusFunction(vGBolus, tGBolus, TGBolus);  % [mmol/min]
    
    % Glucose Infusion
    P.data.GInfusion = zeros(size(P.results.tArray)); % [mmol/min]    

    %% Other
    P = GetCPeptideParameters(P);    
    
    %% Save
    patientSet{ii} = P;
    clear P
end
end

