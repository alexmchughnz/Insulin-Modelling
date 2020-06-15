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
inSimTime = (P.simTime(1) <= P.G{3}.time) & (P.G{3}.time < P.simTime(end));
measTime = minutes(P.G{3}.time(inSimTime) - P.simTime(1));
measG = P.G{3}.value(inSimTime);    


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
    
    % Simulate G(t, d2) with resulting SI.
    P = GIModel(P);  % Required for P2.
    P = IDModel(P);  % Required for QDF.
    P = GCModel(P);
    
    % Find average error G(t, d2) to measured data.
    simG = P.results.G(measTime+1);
    GError = abs(simG - measG)./measG;
    GErrorGrid(ii) = mean(GError);
    
%     figure(1)
%     hold on
%     plot(P.G{3}.time, P.G{3}.value, 'r*')
%     time = P.simTime(1) + P.results.tArray/24/60;
%     plot(time, P.results.G)
    
end

%% Solving
isBest = GErrorGrid == min(GErrorGrid);
originalP.d2 = d2Grid(isBest);


end

