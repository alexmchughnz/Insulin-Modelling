function newPatientSet = MakeOGTTLui(patientSet, showPlots)
% Function for loading OGTTLui data.
% INPUTS:
%   patientSet - cell array of patient structs
% OUTPUT:
%   newPatientSet - updated cell array of patient structs

global CONFIG
DEBUGPLOTS = DebugPlots();
C = LoadConstants();

getrow = @(label, n) repmat(string(label), 1, n) + (0:n-1);

source = "OGTTLui";

%% Load Data
filepath = fullfile(CONFIG.DATAPATH, source, "OGTTLuiMaster.xls");

btOpts = detectImportOptions(filepath, 'Sheet', 'Blood Test');
btOpts.DataRange = 'A3:AY15';
btOpts.VariableNamesRange = 'A2:AY2';
btOpts.RowNamesRange = 'A3:A15';
btOpts = setvartype(btOpts, 'double');
btTable = readtable(filepath, btOpts);
nBtMeas = 9;

pocOpts = detectImportOptions(filepath, 'Sheet', 'Point of Care');
pocOpts.DataRange = 'A3:BX15';
pocOpts.VariableNamesRange = 'A2:BX2';
pocOpts.RowNamesRange = 'A3:A15';
pocOpts = setvartype(pocOpts, 'double');
pocTable = readtable(filepath, pocOpts);
nPocMeas = 15;

infoOpts = detectImportOptions(filepath, 'Sheet', 'Subject Info');
infoOpts.DataRange = 'A2:M14';
infoOpts.VariableNamesRange = 'A1:M1';
infoOpts.RowNamesRange = 'A2:A14';
infoOpts = setvartype(infoOpts, 'double');
infoTable = readtable(filepath, infoOpts);


%% Generate Patients
newPatientSet = {};
for ii = 1:length(patientSet)
    baseP = patientSet{ii};
    assert(baseP.patientNum <= 30, "Invalid patient number.")
    
    allCodes = string(infoTable.Properties.RowNames);
    allNums = arrayfun(@(s) sscanf(s, "%d"), allCodes);
    subpatientCodes = allCodes(allNums == baseP.patientNum);
    hasSubpatients = length(subpatientCodes) > 1;
    
    for ll = 1:length(subpatientCodes)
        P = baseP;
        
        code = subpatientCodes{ll};
        if hasSubpatients
            P.patientCode = P.patientCode + code(end);
            P.patientNum = P.patientNum*100 + ll;
        end
        
        
        %% Patient Info
        P.data.age = infoTable{code, "Age"};  % [years]
        P.data.mass = infoTable{code, "Mass"};  % [kg]
        P.data.height = infoTable{code, "Height"};  % [cm]
        P.data.BMI = infoTable{code, "BMI"};  % [kg/m^2]
        
        %% Trial Times
        btTimes = btTable{code, getrow("TP", nBtMeas)};  % Time of measurement [min]
        pocTimes = pocTable{code, getrow("tG", nPocMeas)};  % Time of measurement [min]
        
        allTimes = [btTimes pocTimes];
        allTimes = allTimes(~isnan(allTimes));
        
        P.data.simTime = [min(allTimes) max(allTimes)];
        P.data.simDuration =  @() floor(diff(P.data.simTime));
        P.results.tArray = (P.data.simTime(1) : P.data.simTime(end))';
        
        
        %% Assay Data
        % Glucose Assay
        isVenous = logical(rmmissing(pocTable{code, getrow('v', nPocMeas)}));
        GPOC = rmmissing(pocTable{code, getrow('G', nPocMeas)});
        tPOC = rmmissing(pocTable{code, getrow('tG', nPocMeas)});
        
        vPOCV = GPOC(isVenous);  % Venous + Test Strip
        tPOCV = tPOC(isVenous);  % [min]
        
        vPOCF = GPOC(~isVenous);  % Finger Prick + Test Strip
        tPOCF = tPOC(~isVenous);  % [min]
        
        vGBT = btTable{code, getrow('G', nBtMeas)};  % Blood Test
        tBT = btTable{code, getrow('TP', nBtMeas)};
        
        % Small workaround for some missing BT data for some subjects.
        isValid = ~isnan(vGBT);
        vGBT = vGBT(isValid);  % [mmol/L]
        tGBT = tBT(isValid);  % [min]
        
        if any(isValid)
            P.data.G.value = vGBT';
            P.data.G.time = tGBT';
        else
            P.data.G.value = GPOC';
            P.data.G.time = tPOC';
        end
        P.data.GFast = @(~) P.data.G.value(1);  % [mmol/L]
        
        % Insulin Assay
        vI = C.pmol2mU(btTable{code, getrow('I', nBtMeas)})';  % [mU/L]
        tI = tBT';  % [min]   
        isValid = ~isnan(vI); 
        
        P.data.I.value = vI(isValid);
        P.data.I.time = tI(isValid);
        
        % C-peptide Assay
        P.data.CPep.value = btTable{code, getrow('C', nBtMeas)}';  % [pmol/L]
        P.data.CPep.time = tBT';  % [min]
        
        %% Trial Inputs
        % Insulin Bolus
        P.data.IType = "human";
        P.data.IDelivery = "subcutaneous";
        
        nIMeas = 2;
        P.data.vIBolus = rmmissing(pocTable{code, getrow('I', nIMeas)});  % [mU]
        P.data.tIBolus = rmmissing(pocTable{code, getrow('tI', nIMeas)});  %  [min]
        P.data.TIBolus = 1;  % [min]
        
        % Glucose Bolus
        P.data.GDelivery = "enteral";
        
        vGBolus = 35;  % [g]
        P.data.vGBolus = vGBolus / C.MGlucose * 1e+3;  % [mmol]
        P.data.tGBolus = 0;  % [min]
        P.data.TGBolus = 1;  % [min]
        
        P = MakeBolusFunctions(P);
        
        % Glucose Infusion
        P.data.GInfusion = zeros(size(P.results.tArray)); % [mmol/min]
        
        
        %% Debug Plots
        if showPlots
            DEBUGPLOTS.MakeOGTTLui = struct();
            DP = DEBUGPLOTS.MakeOGTTLui;
            MakeDebugPlot("OGTTLui Input", P, DP);
            
            plt = plot(tPOCV, vPOCV, 'b*');
            plt.DisplayName = "Venous Test Strip";
            
            plt = plot(tPOCF, vPOCF, 'b+');
            plt.DisplayName = "Finger Prick Test Strip";
            
            if ~isempty(vGBT)
                plt = plot(tGBT, vGBT, 'g*');
                plt.DisplayName = "Blood Test";
            end
            
            xlim([0 inf])
            ylim([4 14])
            
            ylabel("Plasma Glucose [mmol/L]")
            legend()
        end
        
        %% Save
        newPatientSet{end+1} = P;
        clear P
        
    end
end

end
