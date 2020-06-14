function P = SolveSystem(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI


%% Setup
options = odeset('RelTol',1e-5, ...
                 'AbsTol',1e-4);


%% Models
P = GIModel(P, options);
P = IDModel(P, options);
P = GCModel(P, options);

end
