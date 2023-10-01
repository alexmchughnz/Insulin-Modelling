function [P, basisSplines, tExtended, allKnots, knots] = MakeSplineBasisFunctions(P, splineOptions)
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
switch splineOptions.knotType
    case "amount"  % Fixed Number of Knots
        numDataKnots = splineOptions.knots;  
        knots = linspace(tStart, tEnd, numDataKnots)';
        knots = RoundToMultiple(knots, tDelta);
    
    case "location" % Specified Knot Locations
        knots = splineOptions.knots';

    otherwise
        error("Must have valid knotType setting for splines.")
end

% Append extra knots for consistent fitting.
numExtraKnots = splineOptions.order;  % k-th order splines requires k extra knots / 1 spline at each end to fully define all in range.
extraKnotSpacing = round(max(diff(knots)));  % Default width of additional knots.

extraLeftKnots = knots(1) - extraKnotSpacing * [numExtraKnots:-1:1];
extraRightKnots = knots(end) + extraKnotSpacing * [1:+1:numExtraKnots];

allKnots = [extraLeftKnots knots extraRightKnots];

% Splines
tExtended = (allKnots(1) : tDelta : ceil(allKnots(end)))';        % Extended time range covering all splines.

% We require maxSplines of the lowest order splines to interpolate to get numDataKnots fully covered
% by maxOrder splines.
maxOrder = splineOptions.order;
numSplines = @(k) numel(allKnots) - k - 1;    % One spline for each group of k+2 knots.
maxSplines = numSplines(0);
phi = zeros(numel(tExtended), maxSplines, maxOrder+1);  % 3D array of spline functions. Dimensions are [time, spline, order].

%% Spline Creation
% Set up zeroth-order splines.
k = 0;
for ii = 1 : numSplines(k)
    isKnotActive = (allKnots(ii) <= tExtended) & (tExtended < allKnots(ii+1));
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
        prevSplineTerm = (tExtended - allKnots(ii)) / (allKnots(ii+k) - allKnots(ii)) .* phi(:,ii,k);
        nextSplineTerm = (allKnots(ii+k+1) - tExtended) / (allKnots(ii+k+1) - allKnots(ii+1)) .* phi(:,ii+1,k);
        
        phi(:, ii, k+1) = prevSplineTerm + nextSplineTerm;
    end
end

% Return final spline set within range.
basisSplines = phi(:, :, end);  % Of higher order only.

isInTimeRange = (tStart <= tExtended) & (tExtended <= tEnd);  % Crop entries outside data range.
basisSplines = basisSplines(isInTimeRange, :);

isZeroColumn = all(basisSplines == 0, CONST.ROWDIM);  % Delete any columns for support splines not in data range.
basisSplines(:, isZeroColumn) = [];   


% ccSplinePresent = 1 : numSplines(k)-maxOrder; % For kth-order splines, the last k splines are zeroed.


assert(height(basisSplines) == length(P.results.tArray))
% numDataSplines = numDataKnots - maxOrder - 1;  
% assert(width(basisSplines) == numDataKnots)

%% Plottings
plotvars.tSpan = tExtended;
plotvars.phi = phi;
plotvars.knots = allKnots;
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