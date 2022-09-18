function patientSet = MakeDISST(Trial, patientNums)
% Function for loading DISST data.
% INPUTS:
%   Trial       - trial struct
%   patientNums - array of patient nums
% OUTPUT:
%   patientSet - updated cell array of patient structs

CONST = Constants();

%% Load Data
opts = spreadsheetImportOptions(...
    'NumVariables', 9, ...
    'DataRange', 'C3:K53', ...
    'VariableNamesRange', 'C2:K2', ...
    'RowNamesRange', 'B3:B53');
opts = setvartype(opts, 'double');
infoTable = readtable(fullfile(Trial.Config.DATAPATH, Trial.source, "database recent.xls"), opts, ...
    'ReadRowNames', true, ...
    'ReadVariableNames', true);

opts = spreadsheetImportOptions(...
    'NumVariables', 25, ...
    'DataRange', 'A3:Y52', ...
    'VariableNamesRange', '2:2');
opts = setvartype(opts, 'double');
dataTable = readtable(fullfile(Trial.Config.DATAPATH, Trial.source, "DIST recent.xls"), opts, ...
    'ReadRowNames', true, ...
    'ReadVariableNames', true);

nMeas = 5;  % Number of measurements.
    
%% Generate Patients
patientSet = {};
for ii = 1:numel(patientNums)
    P.patientNum = patientNums(ii);
    assert(P.patientNum <= 50, "Invalid patient number.")
    code = dataTable.Properties.RowNames{P.patientNum};
    P.patientCode = code;
    
    %% Patient Info
    P.data.age = infoTable{code, "age_years_"};  % [years]
    P.data.mass = infoTable{code, "weight_kg_"};  % [kg]
    P.data.height = infoTable{code, "height_cm_"};  % [cm]
    P.data.BMI = infoTable{code, "bmi"};  % [kg/m^2]
    
    %% Trial Times
    measTimes = dataTable{code, repmat("time", 1, nMeas) + (1:nMeas)}'/60;  % Time of measurement [min]
    P.data.simTime = [floor(min(measTimes)), ceil(max(measTimes))];
    P.data.simDuration =  floor(diff(P.data.simTime));
    P.results.tArray = (P.data.simTime(1) : 1/60 : P.data.simTime(end))';
    
    %% Assay Data
    % Glucose Assay
    P.data.G.value = dataTable{code, repmat("G", 1, nMeas) + (1:nMeas)}';  % [mmol/L]
    P.data.G.time = measTimes;  % [min]
    P.data.GFast = P.data.G.value(1);  % [mmol/L]
    
    % Insulin Assay
    P.data.I.value = dataTable{code, repmat("I", 1, nMeas) + (1:nMeas)}';  % [mU/L]
    P.data.I.time = measTimes;  % [min]
    
    % C-peptide Assay
    P.data.CPep.value = dataTable{code, repmat("C", 1, nMeas) + (1:nMeas)}';  % [pmol/L]
    P.data.CPep.time = measTimes;  % [min]
    
    %% Trial Inputs
    % Insulin Bolus
    P.data.IType = "human";
    P.data.IDelivery = "intravenous";
    
    P.data.vIBolus = dataTable{code, "IB"} * 1e+3;  % [mU]
    P.data.tIBolus = dataTable{code, "timeIB"}/60;  %  [min]
    P.data.TIBolus = 1;  % [min]
    
    % Glucose Bolus
    P.data.GDelivery = "intravenous";
    
    P.data.vGBolus = dataTable{code, "GB"} / CONST.MGlucose * 1e+3;  % [mmol]
    P.data.tGBolus = dataTable{code, "timeGB"}/60;  % [min]
    P.data.TGBolus = 1;  % [min]
    
    % Glucose Infusion
    P.data.GInfusion = zeros(size(P.results.tArray)); % [mmol/min]    
    
    %% Save
    patientSet{end+1} = P;
    clear P
end
end

