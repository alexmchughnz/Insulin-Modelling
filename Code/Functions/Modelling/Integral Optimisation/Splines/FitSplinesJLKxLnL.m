function [P, A, b, basisSplines] = FitSplinesJLKxLnL(P, numKnots)

CONST = LoadConstants();
GC = P.parameters.GC;

numFixedParameters = 2;

%% Setup
tArray = P.results.tArray;
tMeas = P.data.I.time;
iiMeas = GetTimeIndex(tMeas, tArray);

% Plasma Insulin
[tI, vI] = GetData(P.data.I); % [mU/L]
ppI = griddedInterpolant(tI, vI);  % [mU/L]
I = ppI(tArray);

% Interstital Insulin
Q = GetAnalyticalInterstitialInsulin(I, P);

% Endogenous Secretion
Uen = P.results.Uen;

% Exogenous Insulin
IInput = P.results.IInput;

%% Get Coefficients
% Collect basis functions for splines.
order = 3;
knotLocations = P.data.I.time;

[basisSplines, allKnots] = MakeSplineBasisFunctions(P, order, "knotLocations", knotLocations);
numTotalSplines = size(basisSplines, CONST.ROWWISE);
numTotalParameters = numFixedParameters + numTotalSplines;

% Consider:
% dI/dt = cj*JLK + cx*(1-xL) - cn*nL + kU*Uen - kI*I - kIQ*(I-Q)
%   with nL = sum(nLWeight_i * shape_i).

% Let cWeight_i = shape_i * cn.
% We can express equation as:
% dI/dt = cj*JLK + cx*(1-xL) - sum(cWeight_i * nLWeight_i) - kI*I - kIQ*(I-Q)
cj = IInput/GC.VI;
cx = Uen/GC.VI;
cn = I./(1 + GC.alphaI*I);
cWeights =  basisSplines .* cn;

kI = GC.nK;
kIQ = GC.nI./GC.VI;

%% Integrate I Equation
% I(t) - I(t0) = int{cj}*JLK + int{cx}*(1-xL) - int{sum(cWeight_i * nLWeight_i)} - kI*int{I} - kIQ*int{I-Q}
% Defining CJ = int{cj}, CX = -int{cx}, and CWeights = -int{cWeights}
% CJ*JLK + CX*xL + CWeights*nLWeights = I(t) - I(t0) + kI*int{I} + kIQ*int{I-Q} + CX := C
CJ = cumtrapz(tArray, cj);
CX = -cumtrapz(tArray, cx);
CWeights = -cumtrapz(tArray, cWeights);

intITerm = kI*cumtrapz(tArray, I);
intIQTerm = kIQ*cumtrapz(tArray, I-Q);

I0 = I(1) * ones(size(I));
RHS = [I -I0 intITerm intIQTerm CX];
C = sum(RHS, CONST.ROWWISE);

%% Assemble MLR System
% Extract values at measurement points.
vCJ = CJ(iiMeas);
vCX = CX(iiMeas);
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
% [CJ(t) CX(t) CW_1(t) CW_2(t) ... CW_M(t)] * (JLK; nLW1; nLW2; ... nLWM) = [C(t)] @ measurement times only
CJValues = (vCJ(second) - vCJ(first)) ./ dt;
CXValues = (vCX(second) - vCX(first)) ./ dt;
CWValues = (vCWeights(second,:) - vCWeights(first,:)) ./ dt;
CValues = (vC(second) - vC(first)) ./ dt;

A(:,1) = CJValues;
A(:,2) = CXValues;
A(:,3:numTotalParameters) = CWValues;
b = CValues;

%% Set up Splines
% Enforce absolute value constraints on variables.
lbJLK = 0.01;
ubJLK = 1;
lbxL = 0.6;
ubxL = 0.9;
lbnLW = 0.05;
ubnLW = 0.3;

lb(1) = lbJLK;
lb(2) = lbxL;
lb(3:numTotalParameters) = lbnLW;
ub(1) = ubJLK;
ub(2) = ubxL;
ub(3:numTotalParameters) = ubnLW;

% Enforce constraints on nL spline weights.
numConstraints = numTotalSplines - 1;
[~, vG] = GetData(P.data.G);
GDiffs = vG(2:end) - vG(1:end-1);

for ii = 1:numConstraints
    % Set up difference matrix.
    nLDiffMatrix(ii, ii) = -1;
    nLDiffMatrix(ii, ii+1) = 1;
end

% Directional constraint - delta(nL) should always be opposite to delta(G).
% Expand delta(G) pattern by at edges to constrain "extra" splines.
diffPattern = [GDiffs(1); GDiffs];                                 % One extra left-hand spline.
diffPattern(numel(diffPattern)+1 : numConstraints) = GDiffs(end);  % Any number of extra right-hand splines.
% nL should decrease if G is increasing, thus
% sign(G(ii+1) - G(ii)) * (nL(ii+1) - nL(ii)) < 0
nLDirectionA = sign(diffPattern) .* nLDiffMatrix;
nLDirectionb = zeros(numConstraints, 1);

% Change over time constraint - delta(nL) should be less than
% 0.001 min^-2 (Caumo, 2007).
maxnLRateOfChange = 0.001;  % [min^-2]
dataKnots = allKnots(order + (0:numConstraints))'; % Take off extra knots, keep those that define data range.
tDeltas = diff(dataKnots);
maxnLDeltas = maxnLRateOfChange * tDeltas;
% Absolute change in nL must be less than max change, thus
% abs(nL(ii+1) - nL(ii))  < maxnLDeltas(ii)
% sign(nL(ii+1) - nL(ii)) * (nL(ii+1) - nL(ii)) < maxnLDeltas(ii)
% -sign(G(ii+1) - G(ii)) * (nL(ii+1) - nL(ii)) < maxnLDeltas(ii)
nLChangeA = -sign(diffPattern) .* nLDiffMatrix ;
nLChangeb = maxnLDeltas;

% Place constraints on splines.
% "A" structure is [fixedParams, extraSplines, dataSplines, extraSplines].
iiSplines = numFixedParameters + [1:numTotalSplines];

AConstraint = zeros(2*numConstraints, numTotalParameters);
AConstraint(:, iiSplines) = [nLDirectionA; nLChangeA];

bConstraint = [nLDirectionb; nLChangeb];


% Solve using linear solver.
x = lsqlin(A, b, AConstraint, bConstraint, [], [], lb, ub);

JLK = x(1);
P = ApplyInsulinLossFactor(P, JLK);

P.results.xL = x(2);

nLWeights = x(3:end);
P.results.nL = basisSplines * nLWeights;


%% Plotting
plotvars.basisSplines = basisSplines;
plotvars.nLWeights = nLWeights;
MakePlots(P, plotvars);
end


function MakePlots(P, plotvars)
DP = DebugPlots().FitSplines;

%% Splines
if DP.Splines
    MakeDebugPlot("Splines", P, DP);
    
    % Plot nL.
    plt = plot(P.results.tArray, P.results.nL, 'b');
    plt.DisplayName = "nL";
    
    
    % Plot fitted splines.
    plot(P.results.tArray, plotvars.basisSplines .* plotvars.nLWeights', '--', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    
    
    xlabel("Time [min]")
    ylabel("nL [1/min]")
    
    legend
end

%% nL VS Glucose
if DP.nLGlucose
    persistent nLGlucosePlot;
    
    if isempty(nLGlucosePlot)
        nLGlucosePlot = MakeDebugPlot("nL and G Correlation", P, DP);
    else
        figure(nLGlucosePlot)
    end
    
    
    [tG, vG] = GetData(P.data.G); % [mmol/L]
    iiG = GetTimeIndex(tG, P.results.tArray);
    nL = P.results.nL(iiG);
    nLNorm = nL ./ max(nL);
    
    % Scatter
    sct = plot(vG, nLNorm, 'x');
    sct.DisplayName = "P" + P.patientNum;
    
    % Linear Interpolation
    A(:,1) = vG;
    A(:,2) = ones(size(vG));
    b = nLNorm;
    theta = A\b;
    m = theta(1);
    c = theta(2);
    
    GArray = 4:11;
    plt = plot(GArray, m*GArray+c, 'Color', sct.Color);
    plt.HandleVisibility = 'off';
    
    xlabel("Measured G [mmol/L]")
    ylabel("Normalised nL [-]")
    
    legend
    
end

end

