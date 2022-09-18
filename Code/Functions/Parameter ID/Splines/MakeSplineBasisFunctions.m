function [P, basisSplines, knots] = MakeSplineBasisFunctions(P, splineOptions)
% Creates basis spline functions for a time array.
% This function enforces 'numKnots' knots within the range of tArray, and
% generates additional splines for higher orders.

CONST = Constants();

%% Setup
% Time
tDelta = diff(P.results.tArray(1:2));
tStart = P.results.tArray(1);
tEnd = P.results.tArray(end);

% Knots
maxOrder = splineOptions.order;
numExtraKnots = maxOrder;  % k-th order splines requires k extra knots / 1 spline at each end to fully define all in range.

switch splineOptions.knotType
    case "amount"  % Fixed Number of Knots
    numDataKnots = splineOptions.knots;
    
    knotSpacing = RoundToMultiple((tEnd-tStart)/(numDataKnots-1), tDelta);  % Time width between knots.
    
    knotStart = tStart - numExtraKnots*knotSpacing;  % Time location of knots.
    knotEnd = tEnd + numExtraKnots*knotSpacing;
    knots = knotStart : knotSpacing : knotEnd;
    
    case "location" % Specified Knot Locations
    numDataKnots = numel(splineOptions.knots);
    knots = zeros(1, numDataKnots);
    knots(:) = splineOptions.knots;

    extraKnotSpacing = floor(mean(diff(knots)));  % Default width of additional knots.
    
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

    otherwise
        error("Must have valid knotType setting for splines.")
    
end

% Splines
tSpan = (knots(1) : tDelta : knots(end))';        % Extended time range covering all splines.

% We require maxSplines of the lowest order splines to interpolate to get numDataKnots fully covered
% by maxOrder splines.
numSplines = @(k) numel(knots) - k - 1;    % One spline for each group of k+2 knots.
maxSplines = numSplines(0);
phi = zeros(numel(tSpan), maxSplines, maxOrder+1);  % 3D array of spline functions. Dimensions are [time, spline, order].

%% Spline Creation
% Set up zeroth-order splines.
k = 0;
for ii = 1 : numSplines(k)
    isKnotActive = (knots(ii) <= tSpan) & (tSpan < knots(ii+1));
    iiActive = find(isKnotActive);
    
    if ~isempty(iiActive)
        phi(iiActive, ii, k+1) = ones(size(iiActive));
    end
end

% Recursively define higher order splines.
for k = 1:maxOrder
    % Interpolate N splines of order k-1 to get N-1 splines of order k.
    % The ith spline of order k interpolates the ith and (i+1)th splines of order k-1.
    for ii = 1 : numSplines(k)
        prevSplineTerm = (tSpan - knots(ii)) / (knots(ii+k) - knots(ii)) .* phi(:,ii,k);
        nextSplineTerm = (knots(ii+k+1) - tSpan) / (knots(ii+k+1) - knots(ii+1)) .* phi(:,ii+1,k);
        
        phi(:, ii, k+1) = prevSplineTerm + nextSplineTerm;
    end
end

% Return final spline set within range.
basisSplines = phi(:, :, end);  % Of higher order only.

isInTimeRange = (tStart <= tSpan) & (tSpan <= tEnd);  % Crop entries outside data range.
basisSplines = basisSplines(isInTimeRange, :);

isZeroColumn = all(basisSplines == 0, CONST.COLUMNDIR);  % Delete any columns for support splines not in data range.
basisSplines(:, isZeroColumn) = [];   


% ccSplinePresent = 1 : numSplines(k)-maxOrder; % For kth-order splines, the last k splines are zeroed.


assert(height(basisSplines) == length(P.results.tArray))
% numDataSplines = numDataKnots - maxOrder - 1;  
% assert(width(basisSplines) == numDataKnots)

%% Plottings
plotvars.tSpan = tSpan;
plotvars.phi = phi;
plotvars.knots = knots;
plotvars.order = maxOrder;

P = MakePlots(P, plotvars);
end

function P = MakePlots(P, plotvars)
tag = "MakeSplineBasisFunctions";

%% Splines
    P = AddFigure(P, tag, "BasisSplines");
    
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
    
    xlabel("Time [min]")
end