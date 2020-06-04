function results = SolveSystem(variants)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

%% Setup
options = odeset('RelTol',1e-5, ...
                 'AbsTol',1e-4);
 
tArray = (0 : 0.1 : 150)';
results.tArray = tArray;

%% Models
results = GIModel(results, options, variants);
results = INModel(results, options, variants);
results = GCModel(results, options, variants);

end
