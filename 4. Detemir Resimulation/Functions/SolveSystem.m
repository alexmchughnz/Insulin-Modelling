function P = SolveSystem(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI


%% Setup
options = odeset('RelTol',1e-5, ...
                 'AbsTol',1e-4);
 
tArray = (0 : P.simDuration-1)';  % Simulation time array [min]
P.results.tArray = tArray;


%% Models
P.results = GIModel(P.results, options, variants);
P.results = IDModel(P.results, options, variants);
P.results = GCModel(P.results, options, variants);

end
