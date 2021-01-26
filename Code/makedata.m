% Adapted from "PtDataRead.m".


    
elseif contains(dataset, "DISST")
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
