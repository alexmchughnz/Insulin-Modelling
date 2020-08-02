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

end

