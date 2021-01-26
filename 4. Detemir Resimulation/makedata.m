% Adapted from "PtDataRead.m".
% This code is messy and likely redundant in places: it's a layer to take
% raw data and convert into a .mat format used by the rest of the code.
% As such, there's a whole block for each dataset in the Data directory.

function patientSet = makedata(dataset, patientNums, allowPlots)

load config

global C SC
global DEBUGPLOTS
if ~exist('allowPlots', 'var')
    allowPlots = true;
end


if contains(dataset, "Detemir")
    %% Detemir
    source = "Detemir";
    patientSet = cell(size(patientNums));
    
    % STUB: grab data from previous code structs
    for ii = 1:length(patientNums)
        patientNum = patientNums(ii);
        filename = sprintf("sys%d.mat", patientNum);
        load(fullfile(DATAPATH, source, filename));
        data = sys.Data;
        
        clear P
        P.source = "Detemir";
        P.patientNum = data.PtNo;
        P.patientCode = sprintf("P%d", P.patientNum);
        P.data.mass = data.pt_mass;           % Patient mass [kg]
        
        P.data.trialTime     = [sys.trial_start_t, sys.trial_end_t];
        P.data.trialDuration = @() minutes(diff(P.trialTime));
        
        P.data.simTime     =  [sys.sim_start_t, sys.sim_end_t];
        P.data.simDuration =  @() minutes(diff(P.data.simTime));
        
        P.results.tArray = (0 : P.data.simDuration())';
        
        
        P.data.CPep.value = data.Cpep;        % C-peptide reading [pmol/L]
        P.data.CPep.time = data.Cpep_time;    % Time of C-peptide reading [datetime]
        
        P.data.GOther{1}.value = data.bg1;         % Blood glucose reading [mmol/L]
        P.data.GOther{1}.time = data.bg1_time;     % Time of blood glucose reading [datetime]
        P.data.GOther{2}.value = data.bg2;
        P.data.GOther{2}.time = data.bg2_time;
        P.data.G.value = data.bg3;
        P.data.G.time = data.bg3_time;
        
        P.data.ITotal.value = data.PlasmaI;        % Plasma insulin [pmol/L]
        P.data.ITotal.time  = data.PlasmaI_time;
        
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
            preSimTime = minutes(P.data.simTime(1) - startTime); % How long infusion ran before sim [min]
            
            startTime = 0 - preSimTime;         % Start of infusion [min]
            endTime   = duration - preSimTime;  % End of infusion [min]
            
            % Return infusion data.
            iiInfusion = (startTime <= P.results.tArray) & (P.results.tArray < endTime); % 1 if infusion active [logical]
            P.data.GInfusion = iiInfusion .* MAGIC_DEXTROSE_NUMBER/C.MGlucose/60;  % Glucose infusion over sim time [mmol/min]
        end
        
        % Manual patient insulin peak data.
        if P.patientNum == 1
            P.data.tIPeaks = [115 290 730 1525 1760 2110];
        elseif P.patientNum == 3
            P.data.tIPeaks = [125 670 1542 2153];
        elseif P.patientNum == 4
            P.data.tIPeaks = [152 715 1540 2155];
        end
        
        % Manual patient MSE distribution data.
        if P.patientNum == 1
            P.data.stddevMSE = 6000;
        elseif P.patientNum == 3
            P.data.stddevMSE = 8050;
        elseif P.patientNum == 4
            P.data.stddevMSE = 50;
        end
        
        
        % Save patient structs.
        filename = sprintf("patient%d.mat", P.patientNum);
        save(fullfile(DATAPATH, source, filename), '-struct', 'P');
        fprintf('%s: Saved patient data.\n', P.patientCode);
        
    end
    
    % Generate patient data structs.
    loadpatient = @(n) load(fullfile(DATAPATH, source, sprintf("patient%d.mat", n)));
    for pp = 1 : length(patientNums)
        code = patientNums(pp);
        patientSet{pp} = loadpatient(code);
    end
    
elseif contains(dataset, "DISST")
    %% DISST
    source = "DISST";
    
    % Load table.
    opts = spreadsheetImportOptions(...
        'NumVariables', 9, ...
        'DataRange', 'C3:K53', ...
        'VariableNamesRange', 'C2:K2', ...
        'RowNamesRange', 'B3:B53');
    opts = setvartype(opts, 'double');
    TB = readtable(fullfile(DATAPATH, source, "database recent.xls"), opts, ...
        'ReadRowNames', true, ...
        'ReadVariableNames', true);
    
    opts = spreadsheetImportOptions(...
        'NumVariables', 25, ...
        'DataRange', 'A3:Y52', ...
        'VariableNamesRange', '2:2');
    opts = setvartype(opts, 'double');
    TD = readtable(fullfile(DATAPATH, source, "DIST recent.xls"), opts, ...
        'ReadRowNames', true, ...
        'ReadVariableNames', true);
    
    nMeas = 5;  % Number of measurements.
    pp = 1; % Index for patients saved.
    for ii = patientNums
        code = TD.Properties.RowNames{ii};
        P.source = "DISST";
        P.patientCode = code;
        P.patientNum = ii;
        
        % Patient Info
        P.data.age = TB{code, "age_years_"};
        P.data.BMI = TB{code, "bmi"};
        P.data.mass = TB{code, "weight_kg_"};
        P.data.height = TB{code, "height_cm_"};
        
        % Data
        [P.data.k1, P.data.k2, P.data.k3] = SC.k(P);
        
        
        P.data.G.value = TD{code, repmat("G", 1, nMeas) + (1:nMeas)}';     % Plasma glucose [mmol/L]
        P.data.I.value = TD{code, repmat("I", 1, nMeas) + (1:nMeas)}';     % Plasma insulin [mU/L]
        P.data.CPep.value = TD{code, repmat("C", 1, nMeas) + (1:nMeas)}';  % C-peptide readings [pmol/L]
        
        measTimes = TD{code, repmat("time", 1, nMeas) + (1:nMeas)}'/60;  % Time of measurement [min]
        P.data.G.time = measTimes;
        P.data.I.time = measTimes;
        P.data.CPep.time = measTimes;
        
        %  > Bolus
        vIBolus = TD{code, "IB"} * 1e+3;       % Insulin bolus [mU]
        tIBolus = TD{code, "timeIB"}/60;       % Time of bolus delivery [min]
        TIBolus = 1;                           % Period of bolus action [min]
        % Bolus as function of time, value spread over period.
        % Active if time within period.
        P.data.IBolus = @(t) ((tIBolus <= t) & (t < tIBolus+TIBolus)).*vIBolus/TIBolus;  % [mU/min]
        P.data.tIBolus = tIBolus;
        P.data.vIBolus = vIBolus;
        
        vGBolus = TD{code, "GB"};                % Glucose bolus [g]
        vGBolus = vGBolus / C.MGlucose * 1e+3;  % ''            [mmol]
        tGBolus = TD{code, "timeGB"}/60;         % Time of bolus delivery [min]
        TGBolus = 1;                            % Period of bolus action [min]
        % Bolus as function of time, value spread over period.
        % Active if time within period.
        P.data.GBolus = @(t) ((tGBolus <= t) && (t < tGBolus+TGBolus)).*vGBolus/TGBolus;  % [mmol/min]
        P.data.tGBolus = tGBolus;
        P.data.vGBolus = vGBolus;
        
        % Time
        allTimes = [P.data.CPep.time; P.data.G.time; P.data.I.time];
        P.data.simTime = [floor(min(allTimes)), ceil(max(allTimes))];
        P.data.simDuration =  @() floor(diff(P.data.simTime));
        
        P.results.tArray = (P.data.simTime(1) : 1/60 : P.data.simTime(end))';
        
        % Other Fields
        P.data.GFast = @(~) P.data.G.value(1);
        P.data.GInfusion = zeros(size(P.results.tArray)); % By default, no infusion.
        
        %Other Fields
        stddev = 5/100;
        nTrials = 1000;
        try
            load(ResultsPath(sprintf("%s_montecarlo%gx%d.mat", P.patientCode, stddev, nTrials)), ...
                'stddevError')
            P.data.stddevMSE = stddevError;
        catch
            fprintf("No stddevMSE - run AnalyseInsulinVariance!\n")
            P.data.stddevMSE = 0;
        end
        
        % Save patient structs.
        filename = sprintf("patient%s.mat", P.patientCode);
        save(fullfile(DATAPATH, source, filename), '-struct', 'P');
        fprintf('%s: Saved patient data.\n', P.patientCode);
        
        
        % Generate patient data structs.
        loadpatient = @(code) load(fullfile(DATAPATH, source, sprintf("patient%s.mat", code)));
        if ismember(ii, patientNums)
            patientSet{pp} = loadpatient(P.patientCode);
            pp = pp + 1;
        end
        
        clear P
    end
    
    
elseif contains(dataset, "CREBRF")
    %% CREBRF
    source = "CREBRF";
    
    % Load table.
    opts = spreadsheetImportOptions(...
        'NumVariables', 43, ...
        'DataRange', 'C4:AS47', ...
        'VariableNamesRange', 'C3:AS3', ...
        'RowNamesRange', 'A4:A47');
    opts = setvartype(opts, 'double');
    T = readtable(fullfile(DATAPATH, source, "CREBRFImport.xlsx"), opts, ...
        'ReadRowNames', true, ...
        'ReadVariableNames', true);
    
    % Generate array of patient numbers.
    codes = T.Properties.RowNames;
    nums = zeros(size(codes));
    for pp = 1:length(codes)
        code = codes{pp};
        parts = sscanf(code, "%c%d_");
        nums(pp) = parts(2);
    end
    
    for ii = 1:length(patientNums)
        num = patientNums(ii);
        
        if num == 999
            P.patientCode = "K999ZZ";
            loadpatient = @(code) load(fullfile(DATAPATH, source, sprintf("patient%s.mat", code)));
            P = loadpatient(P.patientCode);
        else
            pp = find(nums == num);
            code = codes{pp};
            
            P.source = "CREBRF";
            P.patientCode = strrep(code, '_', '');
            P.patientNum = num;
            
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
            [P.data.k1, P.data.k2, P.data.k3] = SC.k(P);
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
            
            
        end
        
        %Other Fields
        stddev = 5/100;
        nTrials = 1000;
        try
            load(ResultsPath(sprintf("%s_montecarlo%gx%d.mat", P.patientCode, stddev, nTrials)), ...
                'stddevError')
            P.data.stddevMSE = stddevError;
        catch
            fprintf("No stddevMSE - run AnalyseInsulinVariance!\n")
            P.data.stddevMSE = 0;
        end
        
        % Save patient structs.
        filename = sprintf("patient%s.mat", P.patientCode);
        save(fullfile(DATAPATH, source, filename), '-struct', 'P');
        fprintf('%s: Saved patient data.\n', P.patientCode);
        
        % Generate patient data structs.
        loadpatient = @(code) load(fullfile(DATAPATH, source, sprintf("patient%s.mat", code)));
        patientSet{ii} = loadpatient(P.patientCode);
        
        clear P
    end
    
elseif contains(dataset, "OGTT")
    %% OGTT
    source = "OGTTLui";
    patientSet = cell(size(patientNums));
    allPatientCodes = [];
    
    for ii = 1:length(patientNums)
        %% Load data from CSVs.
        patientNum = patientNums(ii);
        patientCode = sprintf("pt%d", patientNum);
        patientFolder = fullfile(DATAPATH, source, patientCode);
        
        % Pull out each unique subpatient letter.
        filenames = ls(fullfile(patientFolder, patientCode+"*"));
        index = strfind(filenames(1,:), '-') - 1;
        subpatientLetters = unique(filenames(:, index));
        
        for ll = 1:length(subpatientLetters)
            subpatientLabel = patientCode + subpatientLetters(ll);
            
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
                P.patientNum = 1000*patientNum + ll;
                P.patientCode = subpatientLabel;
            else
                P.patientNum = patientNum;
                P.patientCode = patientCode;
            end
            
            P.source = source;
            
            % Patient Info
            %TODO: THIS IS DUMMY DATA - GET REAL DATA FROM LUI
            P.data.age = 25;
            P.data.BMI = 23;
            P.data.mass = 85;
            %             P.data.BSA = T{code, "BSA"};
            P.data.height = 185;
            
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
            P = GetCPepRateParameters(P);
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
            
            P.results.nLxLFitBounds = [];
            
            
            % Save patient structs.
            filename = sprintf("patient%s.mat", P.patientCode);
            save(fullfile(DATAPATH, source, filename), '-struct', 'P');
            fprintf('%s: Saved patient data.\n', P.patientCode);
            
            allPatientCodes = [allPatientCodes, P.patientCode];            
        end        
    end
    
    % Generate patient data structs.
    loadpatient = @(code) load(fullfile(DATAPATH, source, sprintf("patient%s.mat", code)));
    for pp = 1 : length(allPatientCodes)
        code = allPatientCodes(pp);
        patientSet{pp} = loadpatient(code);
    end
end


%% --------------------------------------------------

%% Debug Plots
if allowPlots
    DP = DEBUGPLOTS.makedata;
    if DP.GlucoseInput
        for ii = 1:length(patientSet)
            P = patientSet{ii};
            MakeDebugPlot(P, DP);            
            
            GI = zeros(size(P.results.tArray));
            GB = zeros(size(P.results.tArray));
            for tt = 1:length(P.results.tArray)
                time = P.results.tArray(tt);
                GI(tt) = GetGlucoseDelivery(time, P);
                GB(tt) = P.data.GBolus(time);
            end
            
            subplot(3,1,1)
            plot(P.results.tArray, GI)
            title(sprintf("P%s: G Input", P.patientCode))
            ylabel("[mmol/min]")
            
            subplot(3,1,2)
            plot(P.results.tArray, GB)
            title(sprintf("P%s: G Bolus", P.patientCode))
            ylabel("[mmol/min]")
            
            subplot(3,1,3)
            plot(P.results.tArray, P.data.GInfusion)
            title(sprintf("P%s: G Infusion", P.patientCode))
            ylabel("[mmol/min]")
        end
    end
end

end

function P = GetCPepRateParameters(P)

if P.data.BMI > 30
    CHalfLife1 = 4.55;       % Half-life of C-peptide in compartment 1 [min]
    F = 0.78;
else
    CHalfLife1 = 4.95;       % Half-life of C-peptide in compartment 1 [min]
    F = 0.76;
end

CHalfLife2 = 0.14*P.data.age + 29.2;

a = log(2)/CHalfLife1;
b = log(2)/CHalfLife2;


% Rate constants.
k2 = F*(b-a) + a;
k3 = a*b/(k2); % NOTE: Original code had no factor of 1/2; PDD's thesis does.
k1 = a + b - k2 - k3;
    
P.data.k1 = k1;
P.data.k2 = k2;
P.data.k3 = k3;
end


function durs = time2dur(times)
strs = string(times);
for ii = 1:length(strs)
    if strlength(strs(ii)) < 4
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
