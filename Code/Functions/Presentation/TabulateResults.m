function tables = TabulateResults(patientSet)
% Tabulates patient data.
% INPUTS:
%   patientSet - set of patient structs
% OUTPUT:
%   tables - updated tables

DIM3 = 3;


%% Table Assembly
tables = {};
for ii = 1:length(patientSet)
    P = patientSet{ii};
    code = P.patientCode;
    
    %% Main Results
    name = "MainResults";
    T = table;
    T.Properties.Description = name;
    
    T = AddField(T, code, P.results, "nL", @(x) x(1), "nL [1/min]");
    T = AddField(T, code, P.results, "xL", @(x) x(1), "xL [1]");
    T = AddField(T, code, P.results, "SI", @(x) x*1e+3, "SI [*1e-3 L/mU/min]");
    
    T = AddField(T, code, P.results, "insulinMAPE", @(x) 100*x, "insulinMAPE [%]");
    T = AddField(T, code, P.results, "glucoseMAPE", @(x) 100*x, "glucoseMAPE [%]");
    
    tables = AppendTable(tables, T);
    clear T
    
    
    %% Find Optimal Hepatic Clearance
    T = table;
    T.Properties.Description = 'OptimalHepaticClearance';
    
    T = AddField(T, code, P.data, "stddevMSE");
    T = AddField(T, code, P.results, "minGridMSE");
    T = AddField(T, code, P.results, "minimalErrorRegionSize");
    
    T = AddField(T, code, P.results, "optimalnLRange", @diff);
    T = AddField(T, code, P.results, "optimalnLRange", @(x) x(1), "nLRangeLower");
    T = AddField(T, code, P.results, "optimalnLRange", @(x) x(end), "nLRangeUpper");
    T = AddField(T, code, P.results, "optimalxLRange", @diff);
    T = AddField(T, code, P.results, "optimalxLRange", @(x) x(1), "xLRangeLower");
    T = AddField(T, code, P.results, "optimalxLRange", @(x) x(end), "xLRangeUpper");
    
    tables = AppendTable(tables, T);
    clear T
    
    
    %% Match I Input
    name = "MatchIInput";
    if isfield(P.results, name)        
        if ~exist("IScales", "var")
            IScales = [];
        end
        
        IScales = cat(DIM3, IScales, P.results.MatchIInput.IScales);        
        nLnKScales = P.results.MatchIInput.nLnKScales;
        UenScales = P.results.MatchIInput.UenScales;
    end
    
end

%% Post Processing

% Match I Input
name = "MeanMatchIInput";
T = table;
T.Properties.Description = name;

variableNames = "Uen@"+string(100*UenScales)+"%";
rowNames = "nLnK@"+string(100*nLnKScales)+"%";
for nn = 1:length(nLnKScales)
    T{rowNames(nn), :} = IScales(nn, :);
end
T.Properties.VariableNames = variableNames;

tables = AppendTable(tables, T);

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

function tables = AppendTable(tables, T)
if ~isempty(T)    
    if isempty(tables)
        tables = T;
    else    
        tables = {tables T};
    end    
end

end

