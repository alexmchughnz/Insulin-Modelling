function basisSplines = MakeSplineBasisFunctions(numKnots, order, tArray)
% Creates basis spline functions for a time array.
% This function enforces 'numKnots' knots within the range of tArray, and 
% generates additional splines for higher orders.

PLOTALLSPLINES = true;

%% Setup
% Time
tStart = tArray(1);
tEnd = tArray(end);
tDelta = diff(tArray(1:2));

% Knots
knotSpacing = RoundToMultiple((tEnd-tStart)/(numKnots-1), tDelta);  % Time width of knots.
numExtraKnots = order;  % k-th order splines requires k extra knots/ k-1 splines at each end to fully define all in range.

knotStart = tStart - numExtraKnots*knotSpacing;  % Time location of knots.
knotEnd = tStart + (numKnots+numExtraKnots)*knotSpacing;
knots = knotStart : knotSpacing : knotEnd;   

% Splines
tSplineSpan = (knots(1) : tDelta : knots(end))';        % Extended time range covering all splines.
numSplines = length(knots) - 1;            % One spline in between each adjacent pair of knots.
phi = zeros(length(tSplineSpan), numSplines, order+1);  % 3D array of spline functions. Dimensions are [time, spline, order].

%% Spline Creation
% Set up zeroth-order splines.
k = 0;
for ii = 1 : length(knots)-1
    isKnotActive = (knots(ii) <= tSplineSpan) & (tSplineSpan < knots(ii+1));
    iiActive = find(isKnotActive);
    
    if ~isempty(iiActive)
        phi(iiActive, ii, k+1) = ones(size(iiActive));
    end
end

plotsplines(tSplineSpan, tArray, knots, phi, k, PLOTALLSPLINES);

% Recursively define higher order splines.
for k = 1:order
    % A k-th order spline interpolates the (k-1)-th order spline at the
    % knot k knots away.
    for ii = 1 : numSplines - k 
        prevSplineTerm = (tSplineSpan - knots(ii)) / (knots(ii+k) - knots(ii)) .* phi(:,ii,k);
        nextSplineTerm = (knots(ii+k+1) - tSplineSpan) / (knots(ii+k+1) - knots(ii+1)) .* phi(:,ii+1,k);
        
        phi(:, ii, k+1) = prevSplineTerm + nextSplineTerm;
    end
    
    plotsplines(tSplineSpan, tArray, knots, phi, k, PLOTALLSPLINES);
    
end

% Return final spline set within range.
basisSplines = phi(:, :, end);
isInRange = (tStart <= tSplineSpan) & (tSplineSpan <= tEnd);
basisSplines = basisSplines(isInRange, 1:numKnots);
end


function plotsplines(tSplineSpan, tRange, knots, phi, order, enable)
if enable
    spline = phi(:, :, order+1);
    
    figure
    hold on
    
    ylim([0 1])
    
    plot(tSplineSpan, spline);
    
    knotCoords = [knots; knots];
    line(knotCoords, repmat(ylim', 1, length(knotCoords)), ...
        'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)       
    
    line([tRange(1) tRange(end)], [1 1], ...
        'Color', 'b', 'LineWidth', 2)    
    
    title(sprintf("Order = %d", order))
    xlabel("Time")
    ylabel("Basis Func")
end
end