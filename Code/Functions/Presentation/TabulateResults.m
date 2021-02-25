function tables = TabulateResults(patientSet)
% Tabulates patient data.
% INPUTS:
%   patientSet - set of patient structs
% OUTPUT:
%   tables - updated tables

DIM3 = 3;


%% Table Assembly
tables = {table table table};

for ii = 1:length(patientSet)
    tt = 0;
    P = patientSet{ii};
    code = P.patientCode;
    
    %% Main Results
    tt = tt + 1;
    name = "MainResults";
    tables{tt} = table;
    tables{tt}.Properties.Description = name;
    
    tables{tt} = AddField(tables{tt}, code, P.results, "nL", @(x) x(1), "nL [1/min]");
    tables{tt} = AddField(tables{tt}, code, P.results, "xL", @(x) x(1), "xL [1]");
    tables{tt} = AddField(tables{tt}, code, P.results, "SI", @(x) x*1e+3, "SI [*1e-3 L/mU/min]");
    
    tables{tt} = AddField(tables{tt}, code, P.results, "insulinMAPE", @(x) 100*x, "insulinMAPE [%]");
    tables{tt} = AddField(tables{tt}, code, P.results, "glucoseMAPE", @(x) 100*x, "glucoseMAPE [%]"); 
    
    
    %% Find Optimal Hepatic Clearance
    tt = tt + 1;
    tables{tt} = table;
    tables{tt}.Properties.Description = 'OptimalHepaticClearance';
    
    tables{tt} = AddField(tables{tt}, code, P.data, "stddevMSE");
    tables{tt} = AddField(tables{tt}, code, P.results, "minGridMSE");
    tables{tt} = AddField(tables{tt}, code, P.results, "minimalErrorRegionSize");
    
    tables{tt} = AddField(tables{tt}, code, P.results, "optimalnLRange", @diff);
    tables{tt} = AddField(tables{tt}, code, P.results, "optimalnLRange", @(x) x(1), "nLRangeLower");
    tables{tt} = AddField(tables{tt}, code, P.results, "optimalnLRange", @(x) x(end), "nLRangeUpper");
    tables{tt} = AddField(tables{tt}, code, P.results, "optimalxLRange", @diff);
    tables{tt} = AddField(tables{tt}, code, P.results, "optimalxLRange", @(x) x(1), "xLRangeLower");
    tables{tt} = AddField(tables{tt}, code, P.results, "optimalxLRange", @(x) x(end), "xLRangeUpper");
    
    
    %% Match I Input    
    tt = tt + 1;
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
if isfield(P.results, name)
    tables{tt}.Properties.Description = name;
    
    variableNames = "Uen@"+string(100*UenScales)+"%";
    rowNames = "nLnK@"+string(100*nLnKScales)+"%";
    for nn = 1:length(nLnKScales)
        tables{tt}{rowNames(nn), :} = IScales(nn, :);
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
