% Adapted from "PtDataRead.m".

function patientSet = makedata(dataset, patientNums)

load config

global C SC
global DEBUGPLOTS

if dataset == "Detemir"
    patientSet = cell(size(patientNums));
    
    % STUB: grab data from previous code structs
    for ii = 1:length(patientNums)
        patientNum = patientNums(ii);
        filename = sprintf("sys%d.mat", patientNum);
        load(fullfile(DATAPATH, dataset, filename));
        data = sys.Data;
        
        clear P
        P.source = "Detemir";
        P.patientNum = data.PtNo;
        P.patientCode = sprintf("P%d", P.patientNum);
        P.data.mass = data.pt_mass;           % Patient mass [kg]
        
        P.data.trialTime     = [sys.trial_start_t, sys.trial_end_t];
        P.data.trialDuration = @() minutes(diff(P.trialTime));
        
        P.data.simTime     =  [sys.sim_start_t, sys.sim_end_t];
        P.data.simDuration =  @() minutes(diff(P.data.simTime)) + 1;
        
        P.results.tArray = (0 : P.data.simDuration())';
        P.results.tArray = P.results.tArray(1:end-1);
        
        
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
        save(fullfile(DATAPATH, dataset, filename), '-struct', 'P');
        fprintf('%s: Saved patient data.\n', P.patientCode);
        
    end
    
    % Generate patient data structs.
    loadpatient = @(n) load(fullfile(DATAPATH, dataset, sprintf("patient%d.mat", n)));
    for pp = 1 : length(patientNums)
        n = patientNums(pp);
        patientSet{pp} = loadpatient(n);
    end
    
elseif dataset == "DISST"
    
    % Load table.
    opts = spreadsheetImportOptions(...
        'NumVariables', 7, ...
        'DataRange', 'C3:I53', ...
        'VariableNamesRange', 'C2:I2', ...
        'RowNamesRange', 'B3:B53');
    opts = setvartype(opts, 'double');
    TB = readtable(fullfile(DATAPATH, dataset, "database recent.xls"), opts, ...
        'ReadRowNames', true, ...
        'ReadVariableNames', true);
    
    opts = spreadsheetImportOptions(...
        'NumVariables', 25, ...
        'DataRange', 'A3:Y52', ...
        'VariableNamesRange', '2:2');
    opts = setvartype(opts, 'double');
    TD = readtable(fullfile(DATAPATH, dataset, "DIST recent.xls"), opts, ...
        'ReadRowNames', true, ...
        'ReadVariableNames', true);
    
    N = 5;  % Number of measurements.
    pp = 1; % Index for patients saved.
    for ii = patientNums
        code = TD.Properties.RowNames{ii};
        P.source = "DISST";
        P.patientCode = code;
        P.patientNum = ii;
        
        % Patient Info
        P.data.age = TB{code, "age_years_"};
        P.data.BMI = TB{code, "bmi"};
        
        % Data
        [P.data.k1, P.data.k2, P.data.k3] = SC.k(P);
        
        
        P.data.G.value = TD{code, repmat("G", 1, N) + (1:N)}';     % Plasma glucose [mmol/L]
        P.data.I.value = TD{code, repmat("I", 1, N) + (1:N)}';     % Plasma insulin [mU/L]
        P.data.CPep.value = TD{code, repmat("C", 1, N) + (1:N)}';  % C-peptide readings [pmol/L]
        
        times = TD{code, repmat("time", 1, N) + (1:N)}'/60;  % Time of measurement [min]
        P.data.G.time = times;
        P.data.I.time = times;
        P.data.CPep.time = times;
        
        %  > Bolus
        vIBolus = TD{code, "IB"} * 1e+3;       % Insulin bolus [mU]
        tIBolus = TD{code, "timeIB"}/60;       % Time of bolus delivery [min]
        TIBolus = 1;                          % Period of bolus action [min]
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
        
        % Time
        allTimes = [P.data.CPep.time; P.data.G.time; P.data.I.time];
        P.data.simTime = [floor(min(allTimes)), ceil(max(allTimes))];
        P.data.simDuration =  @() floor(diff(P.data.simTime));
        
        P.results.tArray = (P.data.simTime(1) : 1/60 : P.data.simTime(end))';
        P.results.tArray = P.results.tArray(1:end-1);
        
        % Other Fields
        P.data.GFast = @(~) P.data.G.value(1);
        P.data.GInfusion = zeros(size(P.results.tArray)); % By default, no infusion.
        
        if ismember(ii, patientNums)
            stddev = 5/100;
            nTrials = 1000;
            try
                load(ResultsPath(sprintf("%s_montecarlo%gx%d.mat", P.patientCode, stddev, nTrials)), ...
                    'stddevError')
                P.data.stddevMSE = stddevError;
            catch
                fprintf("No stddevMSE - run AnalyseInsulinVariance!\n")
            end
        end
        
        % Save patient structs.
        filename = sprintf("patient%s.mat", P.patientCode);
        save(fullfile(DATAPATH, dataset, filename), '-struct', 'P');
        fprintf('%s: Saved patient data.\n', P.patientCode);
        
        
        % Generate patient data structs.
        loadpatient = @(code) load(fullfile(DATAPATH, dataset, sprintf("patient%s.mat", code)));
        if ismember(ii, patientNums)
            patientSet{pp} = loadpatient(P.patientCode);
            pp = pp + 1;
        end
        
        clear P
    end
    
elseif dataset == "CREBRF"
    % Load table.
    opts = spreadsheetImportOptions(...
        'NumVariables', 43, ...
        'DataRange', 'C4:AS47', ...
        'VariableNamesRange', 'C3:AS3', ...
        'RowNamesRange', 'A4:A47');
    opts = setvartype(opts, 'double');
    T = readtable(fullfile(DATAPATH, dataset, "CREBRFImport.xlsx"), opts, ...
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
        pp = find(nums == num);
        code = codes{pp};
        
        P.source = "CREBRF";
        P.patientCode = code;
        P.patientNum = num;
        
        % Patient Info
        P.data.age = T{code, "Age"};
        P.data.BMI = T{code, "BMI"};
        
        % Time
        times = [0 2 4 6 8 10 30 60];  % Time of measurement [min]
        N = length(times);
        P.data.simTime = [floor(min(times)), ceil(max(times))];
        P.data.simDuration =  @() floor(diff(P.data.simTime));
        P.results.tArray = (P.data.simTime(1) : 1 : P.data.simTime(end)-1)';
        
        P.data.G.time = times';
        P.data.I.time = times';
        P.data.CPep.time = times';
        
        % Data
        [P.data.k1, P.data.k2, P.data.k3] = SC.k(P);
        
        P.data.G.value = T{code, repmat("G", 1, N) + times}';         % Plasma glucose [mmol/L]
        P.data.I.value = T{code, repmat("I", 1, N) + times}';         % Plasma insulin [uU/mL == mU/L]
        P.data.CPep.value = T{code, repmat("C", 1, N) + times}';      % C-peptide readings [ng/mL]
        P.data.CPep.value = P.data.CPep.value * 1e-9 / 1e-3;          % ''                 [g/L]
        P.data.CPep.value = P.data.CPep.value / C.MCPeptide / 1e-12;  % ''                 [pmol/L]
        
        P.data.IBolus = @(~) 0;  % No fast-I bolus here!
        P.data.GBolus = @(~) 0;  % No G bolus either.
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
        
        %Other Fields
        stddev = 5/100;
        nTrials = 1000;
        try
            load(ResultsPath(sprintf("%s_montecarlo%gx%d.mat", P.patientCode, stddev, nTrials)), ...
                'stddevError')
            P.data.stddevMSE = stddevError;
        catch
            fprintf("No stddevMSE - run AnalyseInsulinVariance!\n")
        end
        
        % Save patient structs.
        filename = sprintf("patient%s.mat", P.patientCode);
        save(fullfile(DATAPATH, dataset, filename), '-struct', 'P');
        fprintf('%s: Saved patient data.\n', P.patientCode);
        
        % Generate patient data structs.
        loadpatient = @(code) load(fullfile(DATAPATH, dataset, sprintf("patient%s.mat", code)));
        patientSet{ii} = loadpatient(P.patientCode);
        
        clear P
    end
end

%% --------------------------------------------------

%% Debug Plots
DP = DEBUGPLOTS.makedata;
if DP.GlucoseInput
    loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
    patientSet = {loadpatient(1), loadpatient(3), loadpatient(4)};
    for ii = 1:length(patientSet)
        P = patientSet{ii};
        MakeDebugPlot(P, DP);
        
        for tt = 1:P.data.simDuration()
            G(tt) = GetGlucoseDelivery(tt, P);
        end
        
        subplot(2,1,1)
        plot(P.results.tArray, G)
        title(sprintf("P%d: G Input", P.patientNum))
        ylabel("[mmol/min]")
        
        subplot(2,1,2)
        plot(P.results.tArray, P.data.GInfusion)
        title(sprintf("P%d: G Infusion", P.patientNum))
        ylabel("[mmol/min]")
    end
end

end
