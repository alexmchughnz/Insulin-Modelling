function patientSet = MakeOGTTLui(Trial, patientNums)
% Function for loading OGTTLui data.
% INPUTS:
%   Trial       - trial struct
%   patientNums - array of patient nums
% OUTPUT:
%   patientSet - updated cell array of patient structs

CONST = Constants();

getrow = @(label, n) repmat(string(label), 1, n) + (0:n-1);


%% Load Data
filepath = fullfile(Trial.Config.DATAPATH, Trial.source, "OGTTLuiMaster.xls");

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
patientSet = {};
for ii = 1:numel(patientNums)
    baseP.patientNum = patientNums(ii);
    baseP.patientCode = "P" + string(baseP.patientNum);
    assert(baseP.patientNum <= 30, "Invalid patient number.")
    
    allCodes = string(infoTable.Properties.RowNames);
    allNums = arrayfun(@(s) sscanf(s, "%d"), allCodes);
    subpatientCodes = allCodes(allNums == baseP.patientNum);
    hasSubpatients = numel(subpatientCodes) > 1;
    
    for ll = 1:numel(subpatientCodes)
        P = baseP;
        
        code = subpatientCodes{ll};
        if hasSubpatients
            P.patientCode = P.patientCode + code(end);
            P.patientNum = Trial.Config.PATIENTSUBNUMBER(P.patientNum, ll);
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
        P.data.simDuration = floor(diff(P.data.simTime));
        P.results.tArray = round(P.data.simTime(1) : P.data.simTime(end))';
        tDelta = diff(P.results.tArray(1:2));
        
        
        %% Assay Data
        % Glucose Assay
        isVenous = logical(rmmissing(pocTable{code, getrow('v', nPocMeas)}));
        vGPOC = rmmissing(pocTable{code, getrow('G', nPocMeas)});
        tGPOC = RoundToMultiple(rmmissing(pocTable{code, getrow('tG', nPocMeas)}), tDelta);
        
        vGPOCV = vGPOC(isVenous);  % Venous + Test Strip
        tGPOCV = tGPOC(isVenous);  % [min]
        
        vGPOCF = vGPOC(~isVenous);  % Finger Prick + Test Strip
        tGPOCF = tGPOC(~isVenous);  % [min]
        
        tBT = btTable{code, getrow('TP', nBtMeas)};        
        isBTValid = ~isnan(tBT);
        tBT = RoundToMultiple(tBT(isBTValid), tDelta);
        
        vGBT = btTable{code, getrow('G', nBtMeas)};  % Blood Test
        vGBT = vGBT(isBTValid);
        
        % Incorporate venous POC Glucose values where BT G not available.
        isGBTMeasured = ~isnan(vGBT);
        vGBT = vGBT(isGBTMeasured);  % [mmol/L]
        tGBT = tBT(isGBTMeasured);  % [min]
        
        assayCV = 3/100; % [percent]
        fingerCV = 7/100; % [percent]
        
        if not(all(isGBTMeasured))
            missingtBT = tBT(~isGBTMeasured);
            
            % Each column is a missing BT timepoint; each row is a
            % replacement venous point. min will find the index of the 
            % closest venous point to each missing BT point.
            distanceMatrix = abs(missingtBT - tGPOCV');
            [~, iiReplace] = min(distanceMatrix);
            
            % Include the unique venous points closest to each missing BT.
            iiReplace = unique(iiReplace);
            
            tFinal = [tGBT tGPOCV(iiReplace)];
            
            [tFinal, order] = sort(tFinal);
            
            vFinal = [vGBT vGPOCV(iiReplace)];
            vFinal = vFinal(order);
            
            P.data.G.value = vFinal';
            P.data.G.time = tFinal';
            
            % Assign correct CVs to different points.
            isAssay = zeros(numel(tFinal), 1);
            isAssay(1:numel(tGBT)) = true;
            isAssay = logical(isAssay(order));
            
            P.data.GCV = zeros(numel(tFinal), 1);
            P.data.GCV(isAssay) = assayCV;
            P.data.GCV(~isAssay) = fingerCV;
            
        else
            P.data.G.value = vGBT';
            P.data.G.time = tBT';
            
            P.data.GCV = assayCV;
            
        end
        
        P.data.GFast = P.data.G.value(1);  % [mmol/L]
        
        % Insulin Assay
        vI = CONST.pmol2mU(btTable{code, getrow('I', nBtMeas)})';  % [mU/L]
        vI = vI(isBTValid);
        
        tI = tBT';  % [min]   
        isBTValid = ~isnan(vI); 
        
        P.data.I.value = vI(isBTValid);
        P.data.I.time = tI(isBTValid);
        
        % C-peptide Assay
        vCPep = btTable{code, getrow('C', nBtMeas)}';  % [pmol/L]
        tCPep = tBT';  % [min]
        
        isBTValid = ~isnan(vCPep); 
        
        P.data.CPep.value = vCPep(isBTValid);
        P.data.CPep.time = tCPep(isBTValid);
        
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
        P.data.vGBolus = vGBolus / CONST.MGlucose * 1e+3;  % [mmol]
        P.data.tGBolus = 0;  % [min]
        P.data.TGBolus = 1;  % [min]
        
        % Glucose Infusion
        P.data.GInfusion = zeros(size(P.results.tArray)); % [mmol/min]
        
        
        %% Plots
        if Trial.Config.ENABLELOADERPLOTS
            tag = "MakeOGTTLui";
            P = AddFigure(P, tag, "Input");
            
            % All Glucose Data
            subplot(2, 1, 1)
            hold on
            
            plt = plot(tGPOCV, vGPOCV, 'b*');
            plt.DisplayName = "Venous Test Strip";
            
            plt = plot(tGPOCF, vGPOCF, 'b+');
            plt.DisplayName = "Finger Prick Test Strip";
            
            if ~isempty(vGBT)
                plt = plot(tGBT, vGBT, 'r*');
                plt.DisplayName = "Blood Test";
            end
            ylim([4 14])
            
            for tt = tBT
                line([tt tt], ylim, 'Color', 'k', 'LineStyle', ':', ...
                    'LineWidth', 0.5, 'HandleVisibility', 'off')
            end
            
            legend()
            
            ylabel("Plasma Glucose [mmol/L]")
            
            % Final Glucose Data
            subplot(2, 1, 2)    
            hold on
            
            plt = plot(P.data.G.time, P.data.G.value, 'g*');
            plt.DisplayName = "Final Profile";
            
            ylim([4 14])
            
            for tt = tBT
                line([tt tt], ylim, 'Color', 'k', 'LineStyle', ':', ...
                    'LineWidth', 0.5, 'HandleVisibility', 'off')
            end
            
            legend()
            ylabel("Plasma Glucose [mmol/L]")
            xlabel("Time [min]")
        end
        
        %% Save
        patientSet{end+1} = P;
        
    end
end

end
