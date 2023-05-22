function [P, A, b, basisSplines] = FitSplinesnL(P, splineOptions)

CONST = Constants();
GC = P.parameters.GC;
minnL = 0;
numFixedParameters = 0;

% Default splineOptions.
defaultSplineOptions.knotType = "location";
defaultSplineOptions.knots = P.data.I.time;
defaultSplineOptions.order = 3;
defaultSplineOptions.maxRate = 0.001;

fields = fieldnames(defaultSplineOptions);
for ii = 1:numel(fields)
    field = fields{ii};
    if ~isfield(splineOptions, field)
        splineOptions.(field) = defaultSplineOptions.(field);
    end
end
P.results.splineOptions = splineOptions;



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
numTotalSplines = size(basisSplines, CONST.COLUMNDIM);
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
C = sum(RHS, CONST.COLUMNDIM);

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
% [CW_1(t) CW_2(t) ... CW_S(t)] * (n_1; n_2; ... n_S) = [C(t)] @ measurement times only
CWValues = (vCWeights(second,:) - vCWeights(first,:)) ./ dt;
CValues = (vC(second) - vC(first)) ./ dt;

A(:,1:numTotalParameters) = CWValues;
b = CValues;

%% Constraints
% 1) Directional constraint - delta(nL) should always be opposite to delta(G).
% 2) Change over time constraint - delta(nL)/delta(t) should be less than 0.001 min^-2 (Caumo, 2007).

% Select knots in data range.
dataKnots = allKnots(splineOptions.order+1 : end-splineOptions.order)';
numDataKnotPairs = numel(dataKnots) - 1;

% At any time t, nL = sum{n_i * bs_i(t)} for all S splines, where bs_i is the i-th basis spline function and n_i that spline's weighting.
% Thus, between adjacent knot times, delta(nL) = sum{n_i * delta(bs_i)} for all splines i
%                                              = [delta(bs_1), delta(bs_2), ..., delta(bs_S)] * (nLWeights)
% For all knot times, nLDiffMatrix = @ knots 1 to 2:  [delta(bs_1), delta(bs_2), ..., delta(bs_S)
%                                    @ knots 2 to 3:   delta(bs_1), delta(bs_2), ..., delta(bs_S)
%                               ...  @ knots N-1 to N: delta(bs_1), delta(bs_2), ..., delta(bs_S)]
% Thus, nLDiffMatrix * (nLWeights) gives a column vector where each entry is delta(nL) for a pair of knots.
for rr = 1 : numDataKnotPairs   % Each ROW represents a pair of knots (time).
    for cc = 1 : numTotalSplines  % Each COL represents a contribution from an individual spline (location).
        t1 = dataKnots(rr);
        t2 = dataKnots(rr+1);
        bs_i = basisSplines(:, cc);

        n1 = SearchArray(t1, tArray);  % Index of knot 1's location in time.
        n2 = SearchArray(t2, tArray);  % Index of knot 2's location in time.

        nLDiffMatrix(rr, cc) = bs_i(n2) - bs_i(n1);  % [N-1 x S]
    end
end

% Constraint 1): nL should decrease if G is increasing, thus for all pairs of knots:
%   delta(G) .* delta(nL)                < 0
%   delta(G) .* nLDiffMatrix*(nLWeights) < 0
[tG, vG] = GetData(P.data.G);
ppG = griddedInterpolant(tG, vG);
GDeltas = diff(ppG(dataKnots));  % [1 x N-1]

nLDirectionA = GDeltas .* nLDiffMatrix;
nLDirectionb = zeros(numDataKnotPairs, 1);

% Constraint 2): delta(nL)/delta(t) should be less than max rate R_max, thus for all pairs of knots:
%   abs(delta(nL)) / delta(t)          <  R_max
%   abs( nLDiffMatrix * (nLWeights) )  <  R_max * (delta(t))
% Simultaneously, we are constraining delta(nL) to opposite direction of delta(G).
% Thus, to ensure each change in xL in (nLDiffMatrix * (nLWeights)) is a positive value, we can rewrite:
%   (sgn(delta(G))) .* nLDiffMatrix * (nLWeights)  < R_max * (delta(t))
tDeltas = diff(dataKnots);
nLChangeA = -sign(GDeltas) .* nLDiffMatrix;  % [N-1 x S]
nLChangeb = splineOptions.maxRate .* tDeltas;

%% Place constraints.
% "A" structure is [fixedParams, extraSplines, dataSplines, extraSplines].
iiSplines = numFixedParameters + [1:numTotalSplines];

AConstraint = zeros(2*numDataKnotPairs, numTotalParameters);
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

