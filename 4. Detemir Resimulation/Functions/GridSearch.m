function originalP = GridSearch(originalP)
% INPUTS:
%   originalP   - patient struct
% OUTPUT:
%   originalP   - modified patient struct with best found d2
    
    
%% Setup
% Input grid.
halfLifeGrid = 10 : 10 : 90;
d2Grid = log(2)./halfLifeGrid;
N = length(d2Grid);

% Results grids.
SIGrid = zeros(originalP.simDuration(), N);
GErrorGrid = zeros(1, N); % Average relative error for each d2 value trialled.

% Measured G (for error comparison)
P = originalP;
[t, G] = GetSimTime(P, P.data.G{3});
inSimTime = (0 <= t) & (t <= P.simDuration()); % [logical]
measTime = t(inSimTime);
measG = G(inSimTime);

%% Search
for ii = 1:N 
    P = originalP;
    fprintf('P%d: Forward simulating with d2=%.4f.\n', ...
            P.patientNum, d2Grid(ii));
    
    % Retrieve d2 value to simulate.
    P.d2 = d2Grid(ii);  % [1/min]
    
    % Fit SI at this d2 value.
    P = FitInsulinSensitivity(P);
    SIGrid(:, ii) = P.SI;
    
    % Simulate G(t, d2) with resulting SI.
    P = GIModel(P);  % Required for P2.
    P = IDModel(P);  % Required for QDF.
    P = GCModel(P);
    
    % Find average error G(t, d2) to measured data.
    simG = P.results.G(measTime+1);
    GError = abs(simG - measG)./measG;
    GErrorGrid(ii) = mean(GError);    
end

%% Solving
isBest = GErrorGrid == min(GErrorGrid);
originalP.d2 = d2Grid(isBest);  % [1/min]


end

