function newPatientSet = MakeOGTTLui(patientSet)
% Function for loading OGTTLui data.
% INPUTS:
%   patientSet - cell array of patient structs
% OUTPUT:
%   newPatientSet - updated cell array of patient structs

global CONFIG DEBUGPLOTS
C = LoadConstants();

source = "OGTTLui";

newPatientSet = {};
for ii = 1:length(patientSet)
    baseP = patientSet{ii};
    
    code = sprintf("pt%d", baseP.patientNum);
    patientFolder = fullfile(CONFIG.DATAPATH, source, code);
    
    % Pull out each unique subpatient letter.
    filenames = ls(fullfile(patientFolder, code+"*"));
    assert(numel(filenames) > 0, "Invalid patient number.");    
    index = strfind(filenames(1,:), '-') - 1;
    subpatientLetters = unique(filenames(:, index));
    
    %% Patient Info
    %TODO: THIS IS DUMMY DATA - GET REAL DATA FROM LUI
    baseP.data.age = 25;  % [years]
    baseP.data.mass = 90;  % [kg]
    baseP.data.height = 185;  % [cm]
    baseP.data.BMI = 23;  % [kg/m^2]
    
    for ll = 1:length(subpatientLetters)
        P = baseP;   
        
        %% Load Data
        subpatientLabel = code + subpatientLetters(ll);        
        metaFile = fullfile(patientFolder, subpatientLabel+"-meta.csv");
        btFile = fullfile(patientFolder, subpatientLabel+"-bt.csv");
        pocFile = fullfile(patientFolder, subpatientLabel+"-poc.csv");
        
        % Load tables.
        opts = detectImportOptions(metaFile, ...
            'NumHeaderLines', 1);
        metaTable = readtable(metaFile, opts, ...
            'ReadVariableNames', true);
        
        opts = detectImportOptions(btFile, ...
            'NumHeaderLines', 1);
        btTable = readtable(btFile, opts, ...
            'ReadVariableNames', true);
        
        opts = detectImportOptions(pocFile, ...
            'NumHeaderLines', 1);
        pocTable = readtable(pocFile, opts, ...
            'ReadVariableNames', true);
        
        if length(subpatientLetters) > 1
            % Workaround subpatients by adding a suffix "00x" to their num.
            P.patientNum = 1000*P.patientNum + ll;
            P.patientCode = P.patientCode + "00"+string(ll);
        end        
        
        %% Trial Times
        % Times here need to be converted from integers representing time
        % of day to an integer representing minutes passed since trial
        % start (e.g. if a trial starts at 1200h, the integer 1315 should
        % be converted to 75 [minutes]).
        allTimes = Time2Duration([metaTable.time; btTable.time; pocTable.time]);
        startTime = min(allTimes);
        endTime = max(allTimes);
        
        time2mins = @(times) minutes(Time2Duration(times) - startTime);
        metaTable.time = time2mins(metaTable.time);
        btTable.time = time2mins(btTable.time);
        pocTable.time = time2mins(pocTable.time);
        
        P.data.simTime = [0, minutes(endTime-startTime)];  % [min]
        P.data.simDuration =  @() floor(diff(P.data.simTime));  % [min]
        P.results.tArray = (P.data.simTime(1) : P.data.simTime(end))';  % [min]
        
        %% Assay Data
        % Glucose Assay
        isValid = ~isnan(btTable.glucose);        
        vBT = btTable.glucose(isValid);  % [mmol/L]
        tBT = btTable.time(isValid);  % [min]
        vPOC = pocTable.glucose;  % [mmol/L]
        tPOC = pocTable.time;  % [min]     
        
        % If duplicated, times, choose blood test.
        isDuplicate = ismember(tPOC, tBT);
        [tG, order] = sort([tBT; tPOC(~isDuplicate)]);
        vG = [vBT; vPOC(~isDuplicate)];
        
        P.data.G.value = vG(order);
        P.data.G.time = tG;
        
        P.data.GFast = @(~) P.data.G.value(1);  % [mmol/L]
        
        % Insulin Assay    
        P.data.I.value = C.pmol2mU(btTable.insulin);  % [mU/L]
        P.data.I.time = btTable.time;  % [min]
        
        % C-peptide Assay        
        P.data.CPep.value = btTable.cpep;  % [pmol/L]
        P.data.CPep.time = btTable.time;  % [min]
        
        %% Trial Inputs
        % Insulin Bolus
        P.data.IType = "human";
        P.data.IDelivery = "subcutaneous";
        
        valid = ~isnan(metaTable.insulin);
        vIBolus = metaTable.insulin(valid) * 1e+3;  % [mU]
        tIBolus = metaTable.time(valid);  % [min]
        TIBolus = 1;  % [min]
        P.data.IBolus = MakeBolusFunction(vIBolus, tIBolus, TIBolus);  % [mU/min]
        
        % Glucose Bolus
        P.data.GDelivery = "enteral";
        
        vGBolus = 35;  % [g]
        vGBolus = vGBolus / C.MGlucose * 1e+3;  % [mmol]
        tGBolus = 0;  % [min]
        TGBolus = 1;  % [min]
        P.data.GBolus = MakeBolusFunction(vGBolus, tGBolus, TGBolus);  % [mmol/min]
        
        % Glucose Infusion
        P.data.GInfusion = zeros(size(P.results.tArray));
        
        %% Other
%         P = GetCPeptideParameters(P);    
        
        %% Debug Plots
        DEBUGPLOTS.MakeOGTTLui = struct();
        DP = DEBUGPLOTS.MakeOGTTLui;
        MakeDebugPlot("OGTTLui Input", P, DP);
        
        plt = plot(tBT, vBT, 'r*');
        plt.DisplayName = "Blood Test";
        
        plt = plot(tPOC, vPOC, 'b*');
        plt.DisplayName = "Point of Care Strip";
        
        ylabel("Plasma Glucose [mmol/L]")
        legend()
    
        %% Save
        newPatientSet = [newPatientSet P];
        clear P
        
    end
end
end

function durs = Time2Duration(times)
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
