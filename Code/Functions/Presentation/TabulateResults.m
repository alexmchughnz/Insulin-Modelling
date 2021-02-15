function T = TabulateResults(T, P)
% Tabulates patient data.
% INPUTS:
%   T - existing table of results
%   P - patient struct
% OUTPUT:
%   T - updated table

code = P.patientCode;
T{code, "code"} = string(code);

T = AddField(T, code, P.results, "nL", @(x) x(1), "nL [1/min]");
T = AddField(T, code, P.results, "xL", @(x) x(1), "xL [1]");
T = AddField(T, code, P.results, "SI", @(x) x*1e+3, "SI [*1e-3 L/mU/min]");

T = AddField(T, code, P.results, "delta2Norm");
T = AddField(T, code, P.results, "delta2NormnL");
T = AddField(T, code, P.results, "delta2NormxL");

T = AddField(T, code, P.results, "insulinMAPE", @(x) 100*x, "insulinMAPE [%]");
T = AddField(T, code, P.results, "glucoseMAPE", @(x) 100*x, "glucoseMAPE [%]");
T = AddField(T, code, P.data, "stddevMSE");
T = AddField(T, code, P.results, "minGridMSE");
T = AddField(T, code, P.results, "minimalErrorRegionSize");

T = AddField(T, code, P.results, "optimalnLRange", @diff);
T = AddField(T, code, P.results, "optimalnLRange", @(x) x(1), "nLRangeLower");
T = AddField(T, code, P.results, "optimalnLRange", @(x) x(end), "nLRangeUpper");
T = AddField(T, code, P.results, "optimalxLRange", @diff);
T = AddField(T, code, P.results, "optimalxLRange", @(x) x(1), "xLRangeLower");
T = AddField(T, code, P.results, "optimalxLRange", @(x) x(end), "xLRangeUpper");

T = AddField(T, code, P.data, "age");
T = AddField(T, code, P.data, "BMI");
T = AddField(T, code, P.data, "mass");
T = AddField(T, code, P.data, "height");
T = AddField(T, code, P.data, "BSA");

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

