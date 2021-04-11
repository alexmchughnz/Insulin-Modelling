function P = FindGutEmptyingRate(P)
% Performs a grid search to find the most suitable gut emptying rate (d2).
% INPUTS:
%   P      - patient struct
% OUTPUT:
%   P   - modified patient struct with best found d2    

if HasPersistent(P, "optimald2")
    P.results.d2 = P.persistents.optimald2;
    return
end
    
%% Setup
% Input grid.
halfLifeGrid = 10 : 10 : 90;
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
    copyP = FitInsulinSensitivity(copyP, false);
    SIGrid(:, ii) = copyP.results.SI;
    
    % Simulate G(t, d2) with resulting SI.
    copyP = SolveSystem(copyP);  % Required for P2 and QDF.
    
    % Find average error G(t, d2) to measured data.
    iiG = GetTimeIndex(tG, P.results.tArray);
    simG = copyP.results.G(iiG);
    
    error = simG - vG;
    error = error(tG >= 0);  % Only evaluate error at true points.
    GErrorGrid(ii) = mean(error.^2); 
end

%% Solving
isBest = GErrorGrid == min(GErrorGrid);
P.results.d2 = d2Grid(isBest);  % [1/min]
P.persistents.optimald2 = P.results.d2;


end

