function [tSpan, splines] = MakeSplineBasisFunctions(nSplines, maxOrder, maxT)


knotWidth = maxT/(nSplines-maxOrder);  % Time width of knots.

knotStart = -maxOrder * knotWidth;
knotEnd = (nSplines + maxOrder - 1) * knotWidth;
knots = knotStart : knotWidth : knotEnd;   % Time location of knots.

tDelta = 0.01;
tSpan = (knots(1) : tDelta : knots(end))';
phi = zeros(length(tSpan), nSplines, maxOrder+1);  % 3D basis func array.

% Set up zeroth-order splines.
order = 0;
for ii = 1 : length(knots)-1
    isKnotActive = (knots(ii) <= tSpan) & (tSpan < knots(ii+1));
    iiActive = find(isKnotActive);
    
    if ~isempty(iiActive)
        phi(iiActive, ii, order+1) = ones(size(iiActive));
    end
end

plotsplines(tSpan, knots, phi, order);

% Recursively define higher order splines.
for order = 1:maxOrder
    for ii = 1 : nSplines + (maxOrder-order)
        prevSplineTerm = (tSpan - knots(ii)) / (knots(ii+order) - knots(ii)) .* phi(:,ii,order);
        nextSplineTerm = (knots(ii+order+1) - tSpan) / (knots(ii+order+1) - knots(ii+1)) .* phi(:,ii+1,order);
        
        phi(:, ii, order+1) = prevSplineTerm + nextSplineTerm;
    end
    
    plotsplines(tSpan, knots, phi, order);
    
end
end


function plotsplines(tSpan, knots, phi, order)
spline = phi(:, :, order+1);

figure(order+1)
hold on

ylim([0 1])

plot(tSpan, spline);

knotCoords = [knots; knots];
line(knotCoords, repmat(ylim', 1, length(knotCoords)), ...
    'Color', 'k', 'LineStyle', ':', 'LineWidth', 0.5)


title(sprintf("Order = %d", order))
xlabel("Time")
ylabel("Basis Func")
end