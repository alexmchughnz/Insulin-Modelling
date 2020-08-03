% Adapted from "PtDataRead.m".

function patientSet = makedata(dataset, patientNums)

load config

global C GC
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
        'NumVariables', 25, ...
        'DataRange', 'A3:Y52', ...
        'VariableNamesRange', '2:2');
    opts = setvartype(opts, 'double');
    
    T = readtable(fullfile(DATAPATH, dataset, "DIST recent.xls"), opts, ...
        'ReadRowNames', true, ...
        'ReadVariableNames', true);
    
    N = 5;  % Number of measurements.
    pp = 1; % Index for patients saved.
    for ii = patientNums
        code = T.Properties.RowNames{ii};
        P.source = "DISST";
        P.patientCode = code;
        P.patientNum = ii;
        
        % Data
        P.data.G.value = T{code, repmat("G", 1, N) + (1:N)}';             % Plasma glucose [mmol/L]
        P.data.I.value = C.mU2pmol(T{code, repmat("I", 1, N) + (1:N)}');  % Plasma insulin [mU/L] -> [pmol/L]
        P.data.CPep.value = T{code, repmat("C", 1, N) + (1:N)}';          % C-peptide readings [pmol/L]
        
        times = T{code, repmat("time", 1, N) + (1:N)}'/60;  % Time of measurement [min]
        P.data.G.time = times;
        P.data.I.time = times;
        P.data.CPep.time = times;
        
        %  > Bolus        
        vIBolus = T{code, "IB"} * 1e+3;       % Insulin bolus [mU]
        tIBolus = T{code, "timeIB"}/60;       % Time of bolus delivery [min]
        TIBolus = 1;                          % Period of bolus action [min]
        % Bolus as function of time, value spread over period.
        % Active if time within period.
        P.data.IBolus = @(t) ((tIBolus <= t) & (t < tIBolus+TIBolus)).*vIBolus/TIBolus;  % [mU/min]
        
        vGBolus = T{code, "GB"};                % Glucose bolus [g]
        vGBolus = vGBolus / C.MGlucose * 1e+3;  % ''            [mmol]
        tGBolus = T{code, "timeGB"}/60;         % Time of bolus delivery [min]
        TGBolus = 1;                            % Period of bolus action [min]
        % Bolus as function of time, value spread over period.
        % Active if time within period.
        P.data.GBolus = @(t) ((tGBolus <= t) && (t < tGBolus+TGBolus)).*vGBolus/TGBolus;  % [mmol/min]
        
        %  > Add early steady-state points.
        earlyTime = -5;
        P.data.I.time = [earlyTime; P.data.I.time];
        P.data.I.value = [P.data.I.value(1); P.data.I.value];
        P.data.G.time = [earlyTime; P.data.G.time];
        P.data.G.value = [P.data.G.value(1); P.data.G.value];
        
        % Time
        P.data.simTime = [earlyTime, max(times)+1];
        P.data.simDuration =  @() floor(diff(P.data.simTime));
        
        P.results.tArray = (earlyTime : 1/60 : P.data.simDuration())';
        P.results.tArray = P.results.tArray(1:end-1);      
        
        % Generate minute-wise insulin profile.
        fakeIData = zeros(size(P.results.tArray));
        
        %  > Interpolate pre-bolus 
        isDataPreBolus = (P.data.I.time <= tIBolus);
        ppI = griddedInterpolant(P.data.I.time(isDataPreBolus), P.data.I.value(isDataPreBolus));  % [pmol/L]
        IInterp = ppI(P.results.tArray);
        
        isSimPreBolus = (P.results.tArray < tIBolus);
        fakeIData(isSimPreBolus) = IInterp(isSimPreBolus);
        
        %  > Bolus
        tAfterIBolus = tIBolus;
        fun = @(x, tdata) P.data.I.value(3) + (tdata > tAfterIBolus).*(x(1)*exp(-x(2)*(tdata - tAfterIBolus)));
        x0 = [1e3; 0.1];
        tdata = P.data.I.time(~isDataPreBolus);
        Idata = P.data.I.value(~isDataPreBolus);
        lb = zeros(size(x0));
        ub = C.mU2pmol(vIBolus)/GC.VI;
        x = lsqcurvefit(fun, x0, tdata, Idata, lb, ub);
        
        fakeIData(~isSimPreBolus) = fun(x, P.results.tArray(~isSimPreBolus));        
        
        %  > Shuffle in fake data points.
        fakeG = vGBolus/GC.VG + P.data.G.value(3);
        [P.data.G.time, order] = sort([P.data.G.time; tGBolus+TGBolus]);
        fakeData = [P.data.G.value; fakeG];
        P.data.G.value = fakeData(order);
        
        fakeI = max(fakeIData);
        [P.data.I.time, order] = sort([P.data.I.time; tAfterIBolus]);
        fakeData = [P.data.I.value; fakeI];
        P.data.I.value = fakeData(order);
        
        
        DP = DEBUGPLOTS.makedata;
        if DP.DISSTBolusFit
            if ismember(P.patientNum, patientNums)
                MakeDebugPlot(P, DP);
                hold on
%                 plt = plot(P.results.tArray, fun(x0, P.results.tArray), ':');
%                 plt.DisplayName = "First Guess";
                plt = plot(P.results.tArray, fakeIData, '.');
                plt.DisplayName = "Fake Minute-wise Data";
                plt = plot(P.data.I.time, P.data.I.value, 'r*');
                plt.DisplayName = "Data";
                fakeii = (P.data.I.value==fakeI);
                plt = plot(P.data.I.time(fakeii), P.data.I.value(fakeii), 'y*');
                plt.DisplayName = "False Point";
                legend()
            end
        end
        
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
                pause(2)
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
