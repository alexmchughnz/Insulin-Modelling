function P = SolveSystem(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P        - patient struct
% OUTPUT:
%   P   - modified patient struct with SI


%% Setup
options = odeset('RelTol',1e-5, ...
    'AbsTol',1e-4);


%% Models
P = GIModel(P, options);

if P.source == "Detemir"
    P = IDModel(P, options);
end
if P.source == "OGTTLui"
    P = SCModel(P, options);
end

P = GCModel(P, options);


%% Fit Evaluation
iiData = GetTimeIndex(P.data.I.time, P.results.tArray);
fitI = P.results.I(iiData);
dataI = P.data.I.value;

MAPE = mean(abs(dataI - fitI)./dataI);
P.results.insulinMAPE = MAPE;
end
