function basisSplines = MakeSplineBasisFunctions(P, order, mode, varargin)
% Creates basis spline functions for a time array.
% This function enforces 'numKnots' knots within the range of tArray, and
% generates additional splines for higher orders.


%% Setup
% Time
tDelta = diff(P.results.tArray(1:2));
tStart = P.results.tArray(1), tDelta;
tEnd = P.results.tArray(end), tDelta;

% Knots
numExtraKnots = order;  % k-th order splines requires k extra knots / k-1 splines at each end to fully define all in range.

if mode == "numKnots"
    % Fixed Number of Knots
    numKnots = varargin{1};
    
    knotSpacing = RoundToMultiple((tEnd-tStart)/(numKnots-1), tDelta, -1);  % Time width of knots.
    
    knotStart = tStart - numExtraKnots*knotSpacing;  % Time location of knots.
    knotEnd = tStart + (numKnots+numExtraKnots)*knotSpacing;
    knots = knotStart : knotSpacing : knotEnd;
    
elseif mode == "knotLocations"
    % Specified Knot Locations
    knots = varargin{1}';    
    extraKnotSpacing = floor(mean(diff(knots)));
    
    % Add knots to cover time range if required.
    if tStart < knots(1)
        knots = [knots(1)-extraKnotSpacing, knots];
    end
    if knots(end) < tEnd
        knots = [knots, knots(end)+extraKnotSpacing];
    end
    
    extraLeftKnots = knots(1) - extraKnotSpacing * [numExtraKnots:-1:1];
    extraRightKnots = knots(end) + extraKnotSpacing * [1:+1:numExtraKnots];
    
    knots = [extraLeftKnots knots extraRightKnots];
    
end

% Splines
tSpan = (knots(1) : tDelta : knots(end))';        % Extended time range covering all splines.
numSplines = length(knots) - 1;                   % One spline in between each adjacent pair of knots.
phi = zeros(length(tSpan), numSplines, order+1);  % 3D array of spline functions. Dimensions are [time, spline, order].

%% Spline Creation
% Set up zeroth-order splines.
k = 0;
for ii = 1 : length(knots)-1
    isKnotActive = (knots(ii) <= tSpan) & (tSpan < knots(ii+1));
    iiActive = find(isKnotActive);
    
    if ~isempty(iiActive)
        phi(iiActive, ii, k+1) = ones(size(iiActive));
    end
end

% Recursively define higher order splines.
for k = 1:order
    % The i-th spline of order k interpolates the (i+k)-th spline
    % of order k-1.
    % A spline of order k connects (k+2) knots.
    for ii = 1 : numSplines - k
        prevSplineTerm = (tSpan - knots(ii)) / (knots(ii+k) - knots(ii)) .* phi(:,ii,k);
        nextSplineTerm = (knots(ii+k+1) - tSpan) / (knots(ii+k+1) - knots(ii+1)) .* phi(:,ii+1,k);
        
        phi(:, ii, k+1) = prevSplineTerm + nextSplineTerm;
    end
end

% Return final spline set within range.
basisSplines = phi(:, :, end);
isInRange = (tStart <= tSpan) & (tSpan <= tEnd);
basisSplines = basisSplines(isInRange, :);

%% Plotting
plotvars.tSpan = tSpan;
plotvars.phi = phi;
plotvars.knots = knots;
plotvars.order = order;
MakePlots(P, plotvars);
end

function MakePlots(P, plotvars)
DP = DebugPlots().MakeSplineBasisFunctions;

%% Splines
if DP.Splines
    MakeDebugPlot("Basis Splines", P, DP);
    
    for k = (0 : plotvars.order) + 1
        subplot(plotvars.order+1, 1, k)
        hold on
        
        spline = plotvars.phi(:, :, k);
        
        plot(plotvars.tSpan, spline);
        
        knotCoords = [plotvars.knots; plotvars.knots];
        line(knotCoords, repmat(ylim', 1, length(knotCoords)), ...
            'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)
        
        line([P.results.tArray(1) P.results.tArray(end)], [1 1], ...
            'Color', 'b', 'LineWidth', 2)
        
        ylim([0 1])
        title(sprintf("Order = %d", k-1))
        ylabel("Spline Value")
    end
    
    xlabel("Time")
end
end