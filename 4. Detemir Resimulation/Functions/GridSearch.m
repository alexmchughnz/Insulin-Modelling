function P = GridSearch(originalP)
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen


global C GC SC
    
    
%% Setup
% Input grid.
halflifeGrid = 10 : 10 : 90;
d2Grid = log(2)./halflifeGrid;
N = length(d2Grid);

% Results grids.
SIGrid = zeros(originalP.simDuration, N);
GErrorGrid = zeros(1, N);


%% Search
for ii = 1:N 
    P = originalP;
    fprintf('P%d: Forward simulating with d2=%.4f.\n', ...
            P.patientNum, d2Grid(ii));
    
    % Retrieve d2 value to simulate.
    P.d2 = d2Grid(ii);
    
    % Fit SI at this d2 value.
    P = FitInsulinSensitivity(P);
    SIGrid(:, ii) = P.SI;
    
    % Simulate G(t) with resulting SI.
    P = GIModel(P);  % Required for P2.
    P = IDModel(P);  % Required for QDF.
    P = GCModel(P);
    
    P.G;
    

%% Solving


end

