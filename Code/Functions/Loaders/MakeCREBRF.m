function patientSet = MakeCREBRF(patientSet)
% Function for loading CREBRF data.
% INPUTS:
%   patientSet  - cell array of patient structs
% OUTPUT:
%   patientSet  - updated cell array of patient structs


global CONFIG
global C

source = "CREBRF";

% Load table.
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

for ii = 1:length(patientSet)
    P = patientSet{ii};
    
    pp = find(nums == P.patientNum);
    assert(numel(pp) == 1, "Invalid patient number.")    
   
    code = codes{pp};
    
    % Patient Info
    P.data.age = T{code, "Age"};
    P.data.BMI = T{code, "BMI"};
    P.data.mass = T{code, "Weight"};
    P.data.BSA = T{code, "BSA"};
    P.data.height = T{code, "Height_cm_"};
    
    % Time
    measTimes = [0 2 4 6 8 10 30];  % Time of measurement [min]
    nMeas = length(measTimes);
    allTimes = [measTimes]';   % Add fake times.
    
    P.data.simTime = [floor(min(allTimes)), ceil(max(allTimes))];
    P.data.simDuration =  @() floor(diff(P.data.simTime));
    P.results.tArray = (P.data.simTime(1) : P.data.simTime(end))';
    
    P.data.G.time = allTimes;
    P.data.I.time = allTimes;
    P.data.CPep.time = allTimes;
    
    % Data
    
    P = GetCPeptideParameters(P);
    
    measG = T{code, repmat("G", 1, nMeas) + measTimes}';  % Plasma glucose [mmol/L]
    P.data.G.value = measG([1:end]);
    measI = T{code, repmat("I", 1, nMeas) + measTimes}';  % Plasma insulin [uU/mL == mU/L]
    P.data.I.value = measI([1:end]);
    measC = T{code, repmat("C", 1, nMeas) + measTimes}';  % C-peptide [ng/mL]
    measC = measC * 1e+3 / 1e-3;                          % ''        [pg/L]
    measC = measC / C.MCPeptide;                          % ''        [pmol/L]
    P.data.CPep.value = measC([1:end]);
    
    P.data.IBolus = @(~) 0;  % No fast-I bolus here!
    P.data.tIBolus = 0;
    P.data.vIBolus = 0;
    
    
    vGBolus = min(0.3*P.data.mass, 30);     % Glucose bolus [g]
    vGBolus = vGBolus / C.MGlucose * 1e+3;  % ''            [mmol]
    tGBolus = 0;                            % Time of bolus delivery [min]
    TGBolus = 1;                            % Period of bolus action [min]
    % Bolus as function of time, value spread over period.
    % Active if time within period.
    P.data.GBolus = @(t) ((tGBolus <= t) && (t < tGBolus+TGBolus)).*vGBolus/TGBolus;  % [mmol/min]
    P.data.tGBolus = tGBolus;
    P.data.vGBolus = vGBolus;
    
    
    P.data.GInfusion = zeros(size(P.results.tArray));  % No glucose infusion in this time range.
    P.data.GFast = @(t) P.data.G.value(1); % Assume starting at fasting.
    
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
    
    
    %% Save
    patientSet{ii} = P;
    clear P
end
end
