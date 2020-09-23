function T = TabulateResults(T, P)
% Tabulates patient data.
% INPUTS:
%   T      - existing table of results
%   P      - patient struct


code = P.patientCode;
T{code, "code"} = string(code);

T = AddField(T, code, P.results, "nL", @(x) x(1));
T = AddField(T, code, P.results, "xL", @(x) x(1));

T = AddField(T, code, P.results, "delta2Norm");
T = AddField(T, code, P.results, "delta2NormnL");
T = AddField(T, code, P.results, "delta2NormxL");

T = AddField(T, code, P.data, "stddevMSE");
T = AddField(T, code, P.results, "minGridMSE");
T = AddField(T, code, P.results, "minimalErrorRegionSize");

T = AddField(T, code, P.results, "optimalnLRange", @diff);
T = AddField(T, code, P.results, "optimalxLRange", @diff);

T = AddField(T, code, P.data, "age");
T = AddField(T, code, P.data, "BMI");
T = AddField(T, code, P.data, "mass");
T = AddField(T, code, P.data, "height");
T = AddField(T, code, P.data, "BSA");

end



function T = AddField(T, code, structP, fieldName, func)
    if ~exist('func', 'var')
        func = @(x) x;
    end

    if isfield(structP, fieldName)    
        T{code, fieldName} = func(structP.(fieldName));
    end
end

