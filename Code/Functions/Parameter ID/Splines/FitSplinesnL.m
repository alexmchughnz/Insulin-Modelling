function [P, A, b, basisSplines] = FitSplinesnL(P, splineOptions)

CONST = Constants();
GC = P.parameters.GC;


minnL = 0;
maxnLRateOfChange = 0.001;  % [min^-2]


numFixedParameters = 0;


%% Setup
tArray = P.results.tArray;
tMeas = P.data.I.time;
iiMeas = SearchArray(tMeas, tArray);

% Plasma Insulin
[tI, vI] = GetData(P.data.I); % [mU/L]
ppI = griddedInterpolant(tI, vI);  % [mU/L]
I = ppI(tArray);

% Interstital Insulin
Q = GetAnalyticalInterstitialInsulin(I, P);

% Endogenous Secretion
Uen = P.results.Uen;

% Exogenous Insulin
Uex = P.results.Uex(P);

%% Get Coefficients
% Collect basis functions for splines.
[P, basisSplines, allKnots] = MakeSplineBasisFunctions(P, splineOptions);
numTotalSplines = size(basisSplines, CONST.ROWDIR);
numTotalParameters = numFixedParameters + numTotalSplines;

% Consider:
% dI/dt = k - cn*nL + kU*Uen - kI*I - kIQ*(I-Q)
%   with nL = sum(nLWeight_i * shape_i).

% Let cWeight_i = shape_i * cn.
% We can express equation as:
% dI/dt = k - sum(cWeight_i * nLWeight_i) + kU*Uen - kI*I - kIQ*(I-Q)
cn = I./(1 + GC.alphaI*I);
cWeights =  basisSplines .* cn;

k = Uex/GC.VI;
kU = (1 - P.results.xL)/GC.VI;
kI = GC.nK;
kIQ = GC.nI./GC.VI;

%% Integrate I Equation
% I(t) - I(t0) = int{k} - int{sum(cWeight_i * nLWeight_i)} + kU*int{Uen} - kI*int{I} - kIQ*int{I-Q}
% Defining CWeights = -int{cWeights}
% CWeights.*nLWeights = I(t) - I(t0) - kU*int{Uen} + kI*int{I} + kIQ*int{I-Q} - int{k} := C
CWeights = -cumtrapz(tArray, cWeights);

intkTerm = cumtrapz(tArray, k);
intUTerm = kU*cumtrapz(tArray, Uen);
intITerm = kI*cumtrapz(tArray, I);
intIQTerm = kIQ*cumtrapz(tArray, I-Q);

I0 = I(1) * ones(size(I));
RHS = [I -I0 -intUTerm intITerm intIQTerm -intkTerm];
C = sum(RHS, CONST.ROWDIR);

%% Assemble MLR System
% Extract values at measurement points.
vCWeights = CWeights(iiMeas, :);
vC = C(iiMeas);

% Find time width of each integral wedge.
first = 1:numel(vC) - 1;
second = 2:numel(vC);

t1 = tMeas(first);
t2 = tMeas(second);
dt = t2 - t1;

% Assemble MLR system by evaluating integral between
% sample points, and normalising by integral width (dt):
% [CW_1(t) CW_2(t) ... CW_M(t)] * (nLW1; nLW2; ... nLWM) = [C(t)] @ measurement times only
CWValues = (vCWeights(second,:) - vCWeights(first,:)) ./ dt;
CValues = (vC(second) - vC(first)) ./ dt;

A(:,1:numTotalParameters) = CWValues;
b = CValues;

%% Set up Splines
% Enforce constraints on nL spline weights.
numConstraints = numTotalSplines - 1;
[tG, vG] = GetData(P.data.G);
ppG = griddedInterpolant(tG, vG);

for ii = 1:numConstraints
    % Set up difference matrix.
    nLDiffMatrix(ii, ii) = -1;
    nLDiffMatrix(ii, ii+1) = 1;
end

% Directional constraint - delta(nL) should always be opposite to delta(G).
% Expand delta(G) pattern by at edges to constrain "extra" splines.
dataKnots = allKnots(splineOptions.order + (0:numConstraints))'; % Take off extra knots, keep those that define data range.
GDiffPattern = diff(ppG(dataKnots));

% nL should decrease if G is increasing, thus
% sign(G(ii+1) - G(ii)) * (nL(ii+1) - nL(ii)) < 0
nLDirectionA = sign(GDiffPattern) .* nLDiffMatrix;
nLDirectionb = zeros(numConstraints, 1);

% Change over time constraint - delta(nL) should be less than
% 0.001 min^-2 (Caumo, 2007).
tDeltas = diff(dataKnots);
maxnLDeltas = maxnLRateOfChange * tDeltas;
% Absolute change in nL must be less than max change, thus
% abs(nL(ii+1) - nL(ii))  < maxnLDeltas(ii)
% sign(nL(ii+1) - nL(ii)) * (nL(ii+1) - nL(ii)) < maxnLDeltas(ii)
% -sign(G(ii+1) - G(ii)) * (nL(ii+1) - nL(ii)) < maxnLDeltas(ii)
nLChangeA = -sign(GDiffPattern) .* nLDiffMatrix ;
nLChangeb = maxnLDeltas;

% Place constraints on splines.
% "A" structure is [fixedParams, extraSplines, dataSplines, extraSplines].
iiSplines = numFixedParameters + [1:numTotalSplines];

AConstraint = zeros(2*numConstraints, numTotalParameters);
AConstraint(:, iiSplines) = [nLDirectionA; nLChangeA];

bConstraint = [nLDirectionb; nLChangeb];


% Solve using linear solver.
lb = minnL * ones(1, numTotalParameters);
x = lsqlin(A, b, AConstraint, bConstraint, [], [], lb);
nLWeights = x;
P.results.nL = basisSplines * nLWeights;


%% Plotting
plotvars.basisSplines = basisSplines;
plotvars.nLWeights = nLWeights;

P = MakePlots(P, plotvars);

end


function P = MakePlots(P, plotvars)
tag = "FitSplines";

%% Splines
    P = AddFigure(P, tag, "Splines");
    
    % Plot nL.
    plt = plot(P.results.tArray, P.results.nL, 'b');
    plt.DisplayName = "nL";
    
    % Plot fitted splines.
    plot(P.results.tArray, plotvars.basisSplines .* plotvars.nLWeights', '--', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    
    ylim([0 0.5])
    
    xlabel("Time [min]")
    ylabel("nL [1/min]")
end

