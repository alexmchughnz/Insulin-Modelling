function tables = TabulateResults(patientSet)
% Tabulates patient data.
% INPUTS:
%   patientSet - set of patient structs
% OUTPUT:
%   tables - updated tables

DIM3 = 3;


tables = {};

%% With All Patients
for ii = 1:numel(patientSet)
    tt = 0;
    P = patientSet{ii};
    code = string(P.patientCode);
    
    %% Main Results
    
    name = "MainResults";
    tt = tt + 1;
    if tt > numel(tables)
        tables{tt} = table;
    end
    tables{tt}.Properties.Description = name;
    
    tables{tt} = AddField(tables{tt}, code, P.results, "JLK", @(x) x, "JLK [/1]");
    
    if all(P.results.nL == P.results.nL(1))
        tables{tt} = AddField(tables{tt}, code, P.results, "nL", @(x) x(1), "nL [1/min]");
    else
        tables{tt} = AddField(tables{tt}, code, P.results, "nL", @(x) mean(x), "mean nL [1/min]");
    end
    tables{tt} = AddField(tables{tt}, code, P.results, "xL", @(x) x(1), "xL [1]");
    
    tables{tt} = AddField(tables{tt}, code, P.results, "d2", @(x) x, "d2 [1/min]");
    tables{tt} = AddField(tables{tt}, code, P.data, "GFast", @(x) x, "Gb [mmol/L]");
    
    if isfield(P.parameters, "SC")
        tables{tt} = AddField(tables{tt}, code, P.parameters.SC, "ks3", @(x) x, "ks3 [1/min]");
    end
    
    if isfield(P.results, "nLGlucose")
        tables{tt} = AddField(tables{tt}, code, P.results.nLGlucose, "coeffs", @(x) x(1), "nL/G slope");
        tables{tt} = AddField(tables{tt}, code, P.results.nLGlucose, "Rsq", @(x) x, "nL/G fit R^2");
    end
    
    tables{tt} = AddField(tables{tt}, code, P.results, "SI", @(x) x*1e+4, "SI [*1e-4 L/mU/min]");
    
    if isfield(P.results, 'fits')
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "insulinMAPE", @(x) 100*x, "insulinMAPE [%]");
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "glucoseMAPE", @(x) 100*x, "glucoseMAPE [%]");
        
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "insulinMSE", @(x) x, "insulinMSE [(mU/L)^2]");
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "glucoseMSE", @(x) x, "glucoseMSE [(mmol/L)^2]");
        
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "insulinMSE", @sqrt, "insulinRMSE [mU/L]");
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "glucoseMSE", @sqrt, "glucoseRMSE [mmol/L]");
    end
    
    if isfield(P, 'persistents')
        tables{tt} = AddField(tables{tt}, code, P.persistents, "stddevMSE");
        tables{tt} = AddField(tables{tt}, code, P.persistents, "stddevSSE");
    end
    
    %% Grid Search Parameters
    name = "GridSearch";
    if isfield(P.results, name)
        tt = tt + 1;
        if tt > numel(tables)
            tables{tt} = table;
        end
        tables{tt}.Properties.Description = name;
        
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "minGridMSE");
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "minGridSSE");
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "minimalErrorRegionSize");
        
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "nLRange", @diff);
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "nLRange", @(x) x(1), "nLRangeLower");
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "nLRange", @(x) x(end), "nLRangeUpper");
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "xLRange", @diff);
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "xLRange", @(x) x(1), "xLRangeLower");
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "xLRange", @(x) x(end), "xLRangeUpper");
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "JLKRange", @diff);
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "JLKRange", @(x) x(1), "JLKLower");
        tables{tt} = AddField(tables{tt}, code, P.results.(name), "JLKRange", @(x) x(end), "JLKUpper");
    end
    
    
end

%% Per Each Patient
tt = numel(tables);
for ii = 1:numel(patientSet)
    P = patientSet{ii};
    
    %% Match I Input - Individual
    name = "MatchIInput";
    if isfield(P.results, name)
        
        % Setup
        if ~exist("allIScales", "var")
            allIScales = [];
        end
        
        nLnKScales = P.results.MatchIInput.nLnKScales;
        UenScales = P.results.MatchIInput.UenScales;
        variableNames = "nLnK@"+string(100*nLnKScales)+"%";
        rowNames = "Uen@"+string(100*UenScales)+"%";
        
        % I Scales
        tt = tt + 1;
        tables{tt} = table;
        tables{tt}.Properties.Description = name +" - " + P.patientCode;
        
        IScales = P.results.MatchIInput.IScales;
        allIScales = cat(DIM3, allIScales, IScales);
        
        for rr = 1:numel(rowNames)
            tables{tt}{rowNames(rr), :} = IScales(rr, :);
        end
        tables{tt}.Properties.VariableNames = variableNames;
        
        
        % I Model Errors
        tt = tt + 1;
        tables{tt} = table;
        tables{tt}.Properties.Description = name +" - " + P.patientCode + " - Error";
        
        IErrors = P.results.MatchIInput.IErrors;
        
        for rr = 1:numel(rowNames)
            tables{tt}{rowNames(rr), :} = IErrors(rr, :);
        end
        tables{tt}.Properties.VariableNames = variableNames;
    end
    
    
    
    %% Match nL - Individual
    name = "MatchnL";
    if isfield(P.results, name)
        
        % Setup
        if ~exist("allIScales", "var")
            allIScales = [];
        end
        
        IInputScales = P.results.MatchnL.IInputScales;
        UenScales = P.results.MatchnL.UenScales;
        variableNames = "IInput@"+string(100*IInputScales)+"%";
        rowNames = "Uen@"+string(100*UenScales)+"%";
        
        % I Scales
        tt = tt + 1;
        tables{tt} = table;
        tables{tt}.Properties.Description = name +" - " + P.patientCode;
        
        IScales = P.results.MatchnL.IScales;
        allIScales = cat(DIM3, allIScales, IScales);
        
        for rr = 1:numel(variableNames)
            tables{tt}{rowNames(rr), :} = IScales(rr, :);
        end
        tables{tt}.Properties.VariableNames = variableNames;
        
        
        % I Model Errors
        tt = tt + 1;
        tables{tt} = table;
        tables{tt}.Properties.Description = name +" - " + P.patientCode + " - Error";
        
        IErrors = P.results.MatchnL.IErrors;
        
        for rr = 1:numel(variableNames)  
            tables{tt}{rowNames(rr), :} = IErrors(rr, :);
        end
        tables{tt}.Properties.VariableNames = variableNames;
    end
end

%% Match I Input - Averaged
name = "MatchIInput";
if isfield(P.results, name)
    tt = tt + 1;
    tables{tt} = table;
    tables{tt}.Properties.Description = name + " - Averaged";
    
    meanIScales = mean(allIScales, DIM3);
    for rr = 1:numel(rowNames)
        tables{tt}{rowNames(rr), :} = meanIScales(rr, :);
    end
    tables{tt}.Properties.VariableNames = variableNames;
end



%% Match nL - Averaged
name = "MatchnL";
if isfield(P.results, name)
    tt = tt + 1;
    tables{tt} = table;
    tables{tt}.Properties.Description = name + " - Averaged";
    
    meanIScales = mean(allIScales, DIM3);
    for rr = 1:numel(variableNames)
        tables{tt}{rowNames(rr), :} = meanIScales(rr, :);
    end
    tables{tt}.Properties.VariableNames = variableNames;
end

end


function T = AddField(T, code, PStruct, fieldName, func, newName)

if ~exist('func', 'var')
    func = @(x) x;
end

if isfield(PStruct, fieldName)
    if exist('newName', 'var')
        T{code, newName} = func(PStruct.(fieldName));
    else
        T{code, fieldName} = func(PStruct.(fieldName));
    end
end

end
