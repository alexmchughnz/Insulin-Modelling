function newPatientSet = MakeOGTTLui(patientSet)
% Function for loading OGTTLui data.
% INPUTS:
%   patientSet  - cell array of patient structs
% OUTPUT:
%   newPatientSet  - updated cell array of patient structs


global CONFIG
global C SC

source = "OGTTLui";

newPatientSet = {};

for ii = 1:length(patientSet)
    baseP = patientSet{ii};
    
    %% Load data from CSVs.
    code = sprintf("pt%d", baseP.patientNum);
    patientFolder = fullfile(CONFIG.DATAPATH, source, code);
    
    % Pull out each unique subpatient letter.
    filenames = ls(fullfile(patientFolder, code+"*"));
    assert(numel(filenames) > 0, "Invalid patient number.");
    
    index = strfind(filenames(1,:), '-') - 1;
    subpatientLetters = unique(filenames(:, index));
    
    % Patient Info
    %TODO: THIS IS DUMMY DATA - GET REAL DATA FROM LUI
    baseP.data.age = 25;
    baseP.data.BMI = 23;
    baseP.data.mass = 85;
    %             P.data.BSA = T{code, "BSA"};
    baseP.data.height = 185;
    
    %% Subpatients
    for ll = 1:length(subpatientLetters)
        P = baseP;
        
        subpatientLabel = code + subpatientLetters(ll);
        
        metaFile = fullfile(patientFolder, subpatientLabel+"-meta.csv");
        btFile = fullfile(patientFolder, subpatientLabel+"-bt.csv");
        pocFile = fullfile(patientFolder, subpatientLabel+"-poc.csv");
        
        % Load tables.
        opts = detectImportOptions(metaFile,...
            'NumHeaderLines', 1);
        metaTable = readtable(metaFile, opts, ...
            'ReadVariableNames', true);
        
        opts = detectImportOptions(btFile,...
            'NumHeaderLines', 1);
        btTable = readtable(btFile, opts, ...
            'ReadVariableNames', true);
        
        opts = detectImportOptions(pocFile,...
            'NumHeaderLines', 1);
        pocTable = readtable(pocFile, opts, ...
            'ReadVariableNames', true);
        
        %% Assemble patient data.
        if length(subpatientLetters) > 1
            % Workaround subpatients by adding a suffix "00x" to their num.
            P.patientNum = 1000*P.patientNum + ll;
            P.patientCode = P.patientCode + "00"+string(ll);
        end
        
        
        
        % Time
        % Times here need to be converted from integers representing time
        % of day to an integer representing minutes passed since trial
        % start (e.g. if a trial starts at 1200h, the integer 1315 should
        % be converted to 75 [minutes]).
        allTimes = time2dur([metaTable.time; btTable.time; pocTable.time]);
        startTime = min(allTimes);
        endTime = max(allTimes);
        
        time2mins = @(times) minutes(time2dur(times) - startTime);
        metaTable.time = time2mins(metaTable.time);
        btTable.time = time2mins(btTable.time);
        pocTable.time = time2mins(pocTable.time);
        
        P.data.simTime = [0, minutes(endTime-startTime)];
        P.data.simDuration =  @() floor(diff(P.data.simTime));
        P.results.tArray = (P.data.simTime(1) : P.data.simTime(end))';
        
        P.data.I.time = btTable.time;
        P.data.CPep.time = btTable.time;
        
        % Data
        P = GetCPeptideParameters(P);
        if isnumeric(btTable.glucose(1)) % Small workaround for some missing BT data for some subjects.
            P.data.G.value = btTable.glucose(~isnan(btTable.glucose));
            P.data.G.time = btTable.time(~isnan(btTable.glucose));
        else
            P.data.G.value = pocTable.glucose;
            P.data.G.time = pocTable.time;
        end
        P.data.I.value = btTable.insulin; % Plasma insulin [mU/L]
        P.data.CPep.value = btTable.cpep; % C-peptide [pmol/L]
        
        % Insulin Bolus
        valid = ~isnan(metaTable.insulin);
        vIBolus = metaTable.insulin(valid) * 1e+3;  % Insulin bolus [mU]
        tIBolus = metaTable.time(valid);            % Time of bolus delivery [min]
        TIBolus = 1;                         % Period of bolus action [min]
        % Bolus as function of time, value spread over period.
        % Active if time within period.
        P.data.IBolus = @(t) dot(vIBolus, ((tIBolus <= t) & (t < tIBolus+TIBolus))) / TIBolus;  % [mU/min]
        P.data.tIBolus = tIBolus;
        P.data.vIBolus = vIBolus;
        
        % Glucose Bolus
        vGBolus = 35;     % Glucose bolus [g]
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
        
        %% Save
        newPatientSet = [newPatientSet P];
        clear P
    end
end
end

function durs = time2dur(times)
strs = string(times);
for ii = 1:length(strs)
    if strlength(strs(ii)) < 4  % e.g. "935" -> "0935"
        strs(ii) = "0"+strs(ii);
    end
end
chararray = char(strs);

addCol = repmat(':', numel(times), 1);
chararray = [chararray(:, 1:2), addCol, chararray(:, 3:4)];

for ii = 1:length(times)
    durs(ii, 1) = duration(chararray(ii, :), 'InputFormat', 'hh:mm');
end
end

