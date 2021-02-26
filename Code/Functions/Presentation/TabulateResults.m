function tables = TabulateResults(patientSet)
% Tabulates patient data.
% INPUTS:
%   patientSet - set of patient structs
% OUTPUT:
%   tables - updated tables

DIM3 = 3;


tables = {};

%% With All Patients
for ii = 1:length(patientSet)
    tt = 0;
    P = patientSet{ii};
    code = P.patientCode;
    
    %% Main Results
    
    name = "MainResults";
    tt = tt + 1;
    tables{tt} = table;
    tables{tt}.Properties.Description = name;
    
    tables{tt} = AddField(tables{tt}, code, P.results, "nL", @(x) x(1), "nL [1/min]");
    tables{tt} = AddField(tables{tt}, code, P.results, "xL", @(x) x(1), "xL [1]");
    tables{tt} = AddField(tables{tt}, code, P.results, "SI", @(x) x*1e+3, "SI [*1e-3 L/mU/min]");
    
    if isfield(P.results, 'fits')
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "insulinMAPE", @(x) 100*x, "insulinMAPE [%]");
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "glucoseMAPE", @(x) 100*x, "glucoseMAPE [%]");
        
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "insulinSSE", @(x) x, "insulinSSE [(mU/L)^2]");
        tables{tt} = AddField(tables{tt}, code, P.results.fits, "glucoseSSE", @(x) x, "glucoseSSE [(mmol/L)^2]");
    end
    
    %% Find Optimal Hepatic Clearance
    
    name = "OptimalHepaticClearance";
    if isfield(P.results, name)
        tt = tt + 1;
        tables{tt} = table;
        tables{tt}.Properties.Description = name;
        
        tables{tt} = AddField(tables{tt}, code, P.data, "stddevMSE");
        tables{tt} = AddField(tables{tt}, code, P.results, "minGridMSE");
        tables{tt} = AddField(tables{tt}, code, P.results, "minimalErrorRegionSize");
        
        tables{tt} = AddField(tables{tt}, code, P.results, "optimalnLRange", @diff);
        tables{tt} = AddField(tables{tt}, code, P.results, "optimalnLRange", @(x) x(1), "nLRangeLower");
        tables{tt} = AddField(tables{tt}, code, P.results, "optimalnLRange", @(x) x(end), "nLRangeUpper");
        tables{tt} = AddField(tables{tt}, code, P.results, "optimalxLRange", @diff);
        tables{tt} = AddField(tables{tt}, code, P.results, "optimalxLRange", @(x) x(1), "xLRangeLower");
        tables{tt} = AddField(tables{tt}, code, P.results, "optimalxLRange", @(x) x(end), "xLRangeUpper");
    end
    
    
end

%% Per Each Patient
tt = 0;
for ii = 1:length(patientSet)
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
        
        for vv = 1:length(variableNames)
            tables{tt}{rowNames(vv), :} = IScales(vv, :);
        end
        tables{tt}.Properties.VariableNames = variableNames;
        
        
        % I Model Errors
        tt = tt + 1;
        tables{tt} = table;
        tables{tt}.Properties.Description = name +" - " + P.patientCode + " - Error";
        
        IErrors = P.results.MatchIInput.IErrors;
        
        for vv = 1:length(variableNames)
            tables{tt}{rowNames(vv), :} = IErrors(vv, :);
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
    for vv = 1:length(variableNames)
        tables{tt}{rowNames(vv), :} = meanIScales(vv, :);
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
