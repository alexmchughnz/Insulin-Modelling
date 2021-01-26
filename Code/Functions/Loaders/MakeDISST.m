function patientSet = MakeDISST(patientNums)
% Function for loading DISST data.
% INPUTS:
%   patientNums - array of patients numbers to load
% OUTPUT:
%   patientSet  - cell array of patient structs 

global CONFIG 
global C SC

source = "DISST";
    
    % Load table.
    opts = spreadsheetImportOptions(...
        'NumVariables', 9, ...
        'DataRange', 'C3:K53', ...
        'VariableNamesRange', 'C2:K2', ...
        'RowNamesRange', 'B3:B53');
    opts = setvartype(opts, 'double');
    TB = readtable(fullfile(CONFIG.DATAPATH, source, "database recent.xls"), opts, ...
        'ReadRowNames', true, ...
        'ReadVariableNames', true);
    
    opts = spreadsheetImportOptions(...
        'NumVariables', 25, ...
        'DataRange', 'A3:Y52', ...
        'VariableNamesRange', '2:2');
    opts = setvartype(opts, 'double');
    TD = readtable(fullfile(CONFIG.DATAPATH, source, "DIST recent.xls"), opts, ...
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
        save(fullfile(CONFIG.DATAPATH, source, filename), '-struct', 'P');
        fprintf('%s: Saved patient data.\n', P.patientCode);
        
        
        % Generate patient data structs.
        loadpatient = @(code) load(fullfile(CONFIG.DATAPATH, source, sprintf("patient%s.mat", code)));
        if ismember(ii, patientNums)
            patientSet{pp} = loadpatient(P.patientCode);
            pp = pp + 1;
        end
        
        clear P
    end
    
    %% Save
    patientSet{ii} = P;
    clear P
end

