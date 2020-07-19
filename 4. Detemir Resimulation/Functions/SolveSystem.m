function P = SolveSystem(P, source)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P        - patient struct
%   source   - string of trial type, e.g. "DISST"
% OUTPUT:
%   P   - modified patient struct with SI


%% Setup
options = odeset('RelTol',1e-5, ...
    'AbsTol',1e-4);


%% Models
if source == "Detemir"
    P = GIModel(P, options);
    P = IDModel(P, options);
    P = GCModel(P, source, options);
else
    P = GIModel(P, options);
    P = GCModel(P, source, options);
end

end
