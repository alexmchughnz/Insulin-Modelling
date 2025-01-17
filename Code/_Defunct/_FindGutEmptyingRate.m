function P = FindGutEmptyingRate(P)
% Performs a grid search to find the most suitable gut emptying rate (d2).
% INPUTS:
%   P      - patient struct
% OUTPUT:
%   P   - modified patient struct with best found d2    

    
%% Setup
% Input grid.
halfLifeGrid = 5 : 10 : 95;
d2Grid = log(2)./halfLifeGrid;
N = length(d2Grid);

% Results grids.
SIGrid = zeros(length(P.results.tArray), N);
GErrorGrid = zeros(1, N); % Average relative error for each d2 value trialled.

% Measured G (for error comparison)
[tG, vG] = GetData(P.data.G);

%% Search
for ii = 1:N   
    d2 = d2Grid(ii);
        
    message = "Searching at d2 = " + string(d2) + "...";
    PrintStatusUpdate(P, message);     
    
    % Retrieve d2 value to simulate.
    copyP = P;
    copyP.results.d2 = d2;  % [1/min]

    % Fit SI at this d2 value.
    copyP = FitInsulinSensitivity(copyP);
    SIGrid(:, ii) = copyP.results.SI;
    
    % Simulate G(t, d2) with resulting SI.
    copyP = SolveSystem(copyP);  % Required for P2 and QDF.
    
    % Find average error G(t, d2) to measured data.
    [~, simG] = GetResultsSample(copyP, tG, copyP.results.G);
    
    glucoseSSE = (simG - vG).^2;
    GErrorGrid(ii) = mean(glucoseSSE); 
end

%% Solving
[~, iiBest] = min(GErrorGrid);
P.results.d2 = d2Grid(iiBest);  % [1/min]
P.persistents.optimald2 = P.results.d2;


end

