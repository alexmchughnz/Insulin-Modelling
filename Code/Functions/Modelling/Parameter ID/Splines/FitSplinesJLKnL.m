function [P, A, b, shapes] = FitSplinesJLKnL(P, numSplines)

CONST = LoadConstants();
GC = P.parameters.GC;

numFixedParameters = 1;
numTotalParameters = numFixedParameters + numSplines;

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
order = 2;
shapes = MakeSplineBasisFunctions(numSplines, order, P.results.tArray);

% Consider:
% dI/dt = cj*JLK - cn*nL + kU*Uen - kI*I - kIQ*(I-Q)
%   with nL = sum(nLWeight_i * shape_i).

% Let cWeight_i = shape_i * cn.
% We can express equation as:
% dI/dt = cj*JLK - sum(cWeight_i * nLWeight_i) + kU*Uen - kI*I - kIQ*(I-Q)
cj = IInput/GC.VI;
cn = I./(1 + GC.alphaI*I);
cWeights =  shapes .* cn;

kU = (1 - P.results.xL)/GC.VI;
kI = GC.nK;
kIQ = GC.nI./GC.VI;

%% Integrate I Equation
% I(t) - I(t0) = int{cj}*JLK - int{sum(cWeight_i * nLWeight_i)} + kU*int{Uen} - kI*int{I} - kIQ*int{I-Q}
% Defining CJ = int{cj} and CWeights = -int{cWeights}
% CJ*JLK + CWeights*nLWeights = I(t) - I(t0) - kU*int{Uen} + kI*int{I} + kIQ*int{I-Q} := C
CJ = cumtrapz(tArray, cj);
CWeights = -cumtrapz(tArray, cWeights);

intUTerm = kU*cumtrapz(tArray, Uen);
intITerm = kI*cumtrapz(tArray, I);
intIQTerm = kIQ*cumtrapz(tArray, I-Q);

I0 = I(1) * ones(size(I));
RHS = [I -I0 -intUTerm intITerm intIQTerm];
C = sum(RHS, CONST.ROWDIR);

%% Assemble MLR System
% Extract values at measurement points.
vCJ = CJ(iiMeas);
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
% [CJ(t) CW_1(t) CW_2(t) ... CW_M(t)] * (JLK; nLW1; nLW2; ... nLWM) = [C(t)] @ measurement times only
CJValues = (vCJ(second) - vCJ(first)) ./ dt;
CWValues = (vCWeights(second,:) - vCWeights(first,:)) ./ dt;
CValues = (vC(second) - vC(first)) ./ dt;

A(:,1) = CJValues;
A(:,2:numTotalParameters) = CWValues;
b = CValues;

%% Solve For Splines
% Enforce constraints on variables.
lbJLK = 0.1;
ubJLK = 1;
lbxLW = 0.05;
ubxLW = 0.3;

lb(1) = lbJLK;
lb(2:numTotalParameters) = lbxLW;
ub(1) = ubJLK;
ub(2:numTotalParameters) = ubxLW;

% Solve using linear solver.
x = lsqlin(A, b, [], [], [], [], lb, ub);
JLK = x(1);

P = ApplyInsulinLossFactor(P, JLK);

nLWeights = x(2:end);
P.results.nL = shapes * nLWeights;


%% Plotting
plotvars.shapes = shapes;
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
    plot(P.results.tArray, plotvars.shapes .* plotvars.nLWeights', '--', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
    
    
    xlabel("Time [min]")
    ylabel("nL [1/min]")
    
    legend
end

end

