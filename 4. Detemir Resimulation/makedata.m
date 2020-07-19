% Adapted from "PtDataRead.m".

function patientSet = makedata(dataset, patientNums)

load config

global C
global DEBUGPLOTS

if dataset == "Detemir"
    patientSet = cell(size(patientNums));
    
    % STUB: grab data from previous code structs
    for ii = 1:length(patientNums)
        patientNum = patientNums(ii);
        filename = sprintf("sys%d.mat", patientNum);
        load(filename);
        data = sys.Data;
        
        clear P
        P.patientNum = data.PtNo;
        P.data.mass = data.pt_mass;           % Patient mass [kg]
        
        P.data.trialTime     = [sys.trial_start_t, sys.trial_end_t];
        P.data.trialDuration = @() minutes(diff(P.trialTime));
        
        P.data.simTime     =  [sys.sim_start_t, sys.sim_end_t];
        P.data.simDuration =  @() minutes(diff(P.data.simTime)) + 1;
        
        P.results.tArray = (0 : P.data.simDuration() - 1)';
        
        P.data.CPep.value = data.Cpep;        % C-peptide reading [pmol/L]
        P.data.CPep.time = data.Cpep_time;    % Time of C-peptide reading [datetime]
        
        P.data.G{1}.value = data.bg1;         % Blood glucose reading [mmol/L?]
        P.data.G{1}.time = data.bg1_time;     % Time of blood glucose reading [datetime]
        P.data.G{2}.value = data.bg2;
        P.data.G{2}.time = data.bg2_time;
        P.data.G{3}.value = data.bg3;
        P.data.G{3}.time = data.bg3_time;
        
        P.data.ITotal.value = data.PlasmaI;        % Plasma insulin [?]
        P.data.ITotal.time  = data.PlasmaI_time;
        
        IBolus = sys.SC.Ibolus;  % Insulin bolus [mU]
        tBolus = sys.SC.T;       % Time of bolus delivery [min]
        TBolus = 5;              % Period of bolus action [min]
        % Bolus as function of time, value spread over period.
        % Active if time within period.
        P.data.IBolus = @(t) ((tBolus <= t) && (t < tBolus+TBolus)).*IBolus/TBolus;
        
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
        save(fullfile(DATAPATH, filename), '-struct', 'P');
        fprintf('P%d: Saved patient data.\n', P.patientNum);
        
    end
    
    % Generate patient data structs.
    loadpatient = @(n) load(fullfile(DATAPATH, sprintf("patient%d.mat", n)));
    for pp = 1 : length(patientNums)
        n = patientNums(pp);
        patientSet{pp} = loadpatient(n);
    end
    
elseif dataset == "DISST"
    
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
