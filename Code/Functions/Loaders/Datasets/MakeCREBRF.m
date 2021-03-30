function patientSet = MakeCREBRF(patientSet)
% Function for loading CREBRF data.
% INPUTS:
%   patientSet  - cell array of patient structs
% OUTPUT:
%   patientSet  - updated cell array of patient structs


global CONFIG
CONST = LoadConstants();

source = "CREBRF";

%% Load Data
opts = spreadsheetImportOptions(...
    'NumVariables', 43, ...
    'DataRange', 'C4:AS47', ...
    'VariableNamesRange', 'C3:AS3', ...
    'RowNamesRange', 'A4:A47');
opts = setvartype(opts, 'double');
T = readtable(fullfile(CONFIG.DATAPATH, source, "CREBRFImport.xlsx"), opts, ...
    'ReadRowNames', true, ...
    'ReadVariableNames', true);

% Generate array of valid patient numbers.
codes = T.Properties.RowNames;
nums = zeros(size(codes));
for pp = 1:length(codes)
    parts = sscanf(codes{pp}, "%c%d_");
    nums(pp) = parts(2);
end

%% Generate Patients
for ii = 1:length(patientSet)
    P = patientSet{ii};
    
    pp = find(nums == P.patientNum);
    assert(numel(pp) == 1, "Invalid patient number.")    
   
    code = codes{pp};
    
    %% Patient Info
    P.data.age = T{code, "Age"};
    P.data.mass = T{code, "Weight"};
    P.data.height = T{code, "Height_cm_"};
    P.data.BMI = T{code, "BMI"};
    
    %% Trial Times
    measTimes = [0 2 4 6 8 10 30]';  % [min]
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
    measC = T{code, repmat("C", 1, nMeas) + measTimes'}';  % [ng/mL]
    measC = measC * 1e+3 / 1e-3;                          % [pg/L]
    measC = measC / CONST.MCPeptide;                          % [pmol/L]
    P.data.CPep.value = measC([1:end]);
    P.data.CPep.time = measTimes;  % [min]   
    
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
    
    % Glucose Bolus
    P.data.GDelivery = "intravenous";
    
    vGBolus = min(0.3*P.data.mass, 30);     % [g]
    P.data.vGBolus = vGBolus / C.MGlucose * 1e+3;  % [mmol]
    P.data.tGBolus = 0;                            % [min]
    P.data.TGBolus = 1;                            % [min]
    
    
    P = MakeBolusFunctions(P);
    
    % Glucose Infusion
    P.data.GInfusion = zeros(size(P.results.tArray));
    
    %% Other
    P = GetCPeptideParameters(P);
    
    %% Save
    patientSet{ii} = P;
    clear P
end
end
