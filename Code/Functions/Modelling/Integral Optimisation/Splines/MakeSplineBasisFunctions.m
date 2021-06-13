function [tSpan, splines] = MakeSplineBasisFunctions(nSplines, maxOrder, maxT)


knotWidth = maxT/(nSplines-maxOrder);  % Time width of knots.

knotStart = -maxOrder * knotWidth;
knotEnd = (nSplines + maxOrder - 1) * knotWidth;
knots = knotStart : knotWidth : knotEnd;   % Time location of knots.

tDelta = 0.01;
tSpan = (knots(1) : tDelta : knots(end))';
phi = zeros(length(tSpan), nSplines, maxOrder+1);  % 3D basis func array.

% Set up zeroth-order splines.
order = 1;
for ii = 1 : length(knots)-1
    isKnotActive = (knots(ii) <= tSpan) & (tSpan < knots(ii+1));
    iiActive = find(isKnotActive);
    
    if ~isempty(iiActive)
        phi(iiActive, ii, order) = ones(size(iiActive));
    end
end

splines = phi(:,:,1);
plot(tSpan, splines)