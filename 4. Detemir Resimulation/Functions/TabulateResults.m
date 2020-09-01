function T = TabulateResults(T, P)
% Tabulates patient data.
% INPUTS:
%   T      - existing table of results
%   P      - patient struct

code = P.patientCode;

T{code, "nL"} = P.results.nL(1);
T{code, "xL"} = P.results.xL(1);
if isfield(P.results, "delta2Norm")
    T{code, "delta2Norm"} = P.results.delta2Norm;
end
if isfield(P.results, "delta2NormnL")
    T{code, "delta2NormnL"} = P.results.delta2NormnL;
end
if isfield(P.results, "delta2NormxL")
    T{code, "delta2NormxL"} = P.results.delta2NormxL;
end
if isfield(P.data, "stddevMSE")
    T{code, "stddevMSE"} = P.data.stddevMSE;
    if isfield(P.results, "minGridMSE")
        T{code, "minGridMSE"} = P.results.minGridMSE;
        T{code, "SD min error %"} = 100 * P.data.stddevMSE / P.results.minGridMSE;
    end
end
if isfield(P.results, "minimalErrorRegionSize")    
    T{code, "MERSize"} = P.results.minimalErrorRegionSize;
end
if isfield(P.results, "optimalnLRange")    
    T{code, "optnLRange"} = diff(P.results.optimalnLRange);
    T{code, "optxLRange"} = diff(P.results.optimalxLRange);
end

end

