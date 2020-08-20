function P = FindGutEmptyingRate(P)
% Performs a grid search to find the most suitable gut emptying rate (d2).
% INPUTS:
%   P      - patient struct
% OUTPUT:
%   P   - modified patient struct with best found d2
    
    
%% Setup
% Input grid.
halfLifeGrid = 10 : 10 : 90;
d2Grid = log(2)./halfLifeGrid;
N = length(d2Grid);

% Results grids.
SIGrid = zeros(length(P.results.tArray), N);
GErrorGrid = zeros(1, N); % Average relative error for each d2 value trialled.

% Measured G (for error comparison)
[tG, vG] = GetSimTime(P, P.data.G);

%% Search
fprintf('P%d: Searching for best d2...\n', ...
        P.patientNum);
for ii = 1:N 
    copyP = P;
    
    % Retrieve d2 value to simulate.
    copyP.results.d2 = d2Grid(ii);  % [1/min]
    
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


end

