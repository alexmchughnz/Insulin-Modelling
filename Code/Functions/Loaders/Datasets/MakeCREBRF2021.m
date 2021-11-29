function patientSetOut = MakeCREBRF2021(Trial, patientNums)
% Function for loading CREBRF2021 data.
% INPUTS:
%   T           - trial struct
%   patientNums - patient numbers to load
% OUTPUT:
%   patientSet  - updated cell array of patient structs

CONST = LoadConstants();

%% Load Data
opts = spreadsheetImportOptions(...
    'NumVariables', 20, ...
    'DataRange', 'B9:U197', ...
    'VariableNamesRange', 'B8:U8', ...
    'RowNamesRange', 'A9:A197');
opts = setvartype(opts, 'double');
T = readtable(fullfile(Trial.Config.DATAPATH, Trial.source, "CREBRF2021Import.xlsx"), opts, ...
    'ReadRowNames', true, ...
    'ReadVariableNames', true);

% Generate array of patient numbers in spreadsheet.
loadCodes = string(T.Properties.RowNames);
loadNums = zeros(size(loadCodes));
for ii = 1:length(loadCodes)
    loadNums(ii) = str2double(regexp(loadCodes{ii}, "\d*", "match", "once"));
end

%% Generate Patients
% Get codes for all numbers and initialise patient set.
% If number is duplicated, just load both codes.
patientNums = unique(patientNums);
patientSet = {};
for ii = 1:numel(patientNums)
    num = patientNums(ii);
    patientIndex = find(loadNums == num);
    assert(numel(patientIndex) > 0, "Invalid patient number: " + string(num))
    
    for pp = 1:numel(patientIndex)
        P = struct();
        P.patientNum = num + 1000*(pp-1); % Add an offset for duplicate numbers.
        P.patientCode = loadCodes(patientIndex(pp));
        patientSet{end+1} = P;
    end
end

patientSetOut = {};
for ii = 1:numel(patientSet)
    P = patientSet{ii};
    code = P.patientCode;
    
    %% Patient Info
    P.data.age = T{code, "AGE"};
    P.data.mass = T{code, "WEIGHT"};
    P.data.height = T{code, "HEIGHT"};
    P.data.BMI = T{code, "BMI"};
    P.data.BSA = T{code, "BSA"};
    
    %% Trial Times
    measTimes = [0 30 60 150]';  % [min]
    nMeas = length(measTimes);
    
    P.data.simTime = [floor(min(measTimes)), ceil(max(measTimes))];
    P.data.simDuration =  floor(diff(P.data.simTime));
    P.results.tArray = (P.data.simTime(1) : P.data.simTime(end))';  % [min]
    
    %% Assay Data
    % Glucose Assay
    measG = T{code, repmat("G", 1, nMeas) + measTimes'}';  % [mmol/L]
    P.data.G.value = measG([1:end]);  % [mmol/L]
    P.data.G.time = measTimes;  % [min]
    P.data.GFast = P.data.G.value(1);  % [mmol/L]
    
    % Insulin Assay
    measI = T{code, repmat("I", 1, nMeas) + measTimes'}';  % [uU/mL == mU/L]
    P.data.I.value = measI([1:end]);  % [mU/L]
    P.data.I.time = measTimes;  % [min]
    
    % C-peptide Assay
    measTimesC = measTimes(1:end-1);  % one fewer CPep measurement
    nMeasC = length(measTimesC);
    measC = T{code, repmat("C", 1, nMeasC) + measTimesC'}';  % [ng/mL]
    measC = measC * 1e+3 / 1e-3;             % [pg/L]
    measC = measC / CONST.MCPeptide;         % [pmol/L]
    P.data.CPep.value = measC([1:end]);
    P.data.CPep.time = measTimesC;           % [min]
    
    % Clear NaNs
    GNaN = isnan(P.data.G.value);
    INaN = isnan(P.data.I.value);
    CNaN = isnan(P.data.CPep.value);
    
    P.data.G.value = P.data.G.value(~GNaN);
    P.data.I.value = P.data.I.value(~INaN);
    P.data.CPep.value = P.data.CPep.value(~CNaN);
    P.data.G.time = P.data.G.time(~GNaN);
    P.data.I.time = P.data.I.time(~INaN);
    P.data.CPep.time = P.data.CPep.time(~CNaN);
    
    %% Trial Inputs
    % Insulin Bolus
    P.data.IType = "none";
    P.data.IDelivery = "none";
    P.data.vIBolus = 0;
    P.data.tIBolus = 0;                            % [min]
    P.data.TIBolus = 0;                            % [min]
    
    % Glucose Bolus
    P.data.GDelivery = "intravenous";
    
    vGBolus = 0;     % [g]  PLACEHOLDER - will need to add mixed meal
    P.data.vGBolus = vGBolus / CONST.MGlucose * 1e+3;  % [mmol]
    P.data.tGBolus = 0;                            % [min]
    P.data.TGBolus = 10;                            % [min]
    
    % Glucose Infusion
    P.data.GInfusion = zeros(size(P.results.tArray));
    
    %% Validation
    valid = true;
    
    valid = valid && (numel(P.data.CPep.value) >= 4);
    valid = valid && (numel(P.data.I.value) >= 3);
    
    %% Save
    if valid
        P.patientCode = upper(strrep(P.patientCode, "_", "-"));
        patientSetOut{end+1} = P;
    end
end
end
