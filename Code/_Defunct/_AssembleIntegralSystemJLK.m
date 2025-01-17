function [A, b, IFunc, QFunc, paramNames] = AssembleIntegralSystemJLK(P, I, Q)

CONST = LoadConstants();
GC = P.parameters.GC;

paramNames = ["JLK"];

%% Setup
tArray = P.results.tArray;
tMeas = P.data.I.time;
iiMeas = GetTimeIndex(tMeas, tArray);

% Plasma Insulin
if ~exist("I", "var")
    [tI, vI] = GetData(P.data.I); % [mU/L]
    ppI = griddedInterpolant(tI, vI);  % [mU/L]
    I = ppI(tArray);
end

% Interstital Insulin
if ~exist("Q", "var")    
    Q = GetAnalyticalInterstitialInsulin(I, P);
end

% Endogenous Secretion
Uen = P.results.Uen;

% Exogenous Insulin
IInput = P.results.IInput;

%% Get Coefficients
% Consider dI/dt = cj*JLK - kIaI*(I/(1+aI*I)) + kU*Uen - kI*I - kIQ*(I-Q):
cj = IInput/GC.VI;

kIaI = P.results.nL;
kU = (1 - P.results.xL)/GC.VI;
kI = GC.nK;
kIQ = GC.nI./GC.VI;

% Also consider dQ/dt = -cQ*Q + cI*I:
cQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q - easier to use
cI = GC.nI/GC.VQ;  % Constant term coefficent of I - easier to use

%% Integrate I Equation
% I(t) - I(t0) = int{cj}*JLK - kIaI*int{I/(1+aI*I)} + kU*int{Uen} - kI*int{I} - kIQ*int{I-Q}
% Defining CJ = int{cj}
% CJ*JLK = I(t) - I(t0) + kIaI*int{I/(1+aI*I)} - kU*int{Uen} + kI*int{I} + kIQ*int{I-Q} := C
CJ = cumtrapz(tArray, cj);

intIaITerm = kIaI.*cumtrapz(tArray, I./(1 + GC.alphaI*I));
intUTerm = kU*cumtrapz(tArray, Uen);
intITerm = kI*cumtrapz(tArray, I);
intIQTerm = kIQ*cumtrapz(tArray, I-Q);

I0 = I(1) * ones(size(I));
RHS = [I -I0 intIaITerm -intUTerm intITerm intIQTerm];
C = sum(RHS, CONST.COLUMNDIM);

%% Make Minute-Wise Q and I Functions
% I(t) = I(t0) + CJ*JLK + CN*nL + kU*int{Uen} - kI*int{I} - kIQ*int{I-Q}
IFunc = @(JLK, nL, I, Q) I(1) + CJ*JLK + CN*nL + intUTerm - intITerm - intIQTerm;

% Q(t) = Q(t0) - cQ*int{Q} + cI*int{I}
QFunc = @(I, Q) Q(1) - cQ*cumtrapz(tArray, Q) + cI*cumtrapz(tArray, I);

%% Assemble MLR System
% Extract values at measurement points.
vCJ = CJ(iiMeas);
vC = C(iiMeas);

% Find time width of each integral wedge. 
first = 1:numel(vC) - 1;
second = 2:numel(vC);

t1 = tMeas(first);
t2 = tMeas(second);
dt = t2 - t1; 

% Assemble MLR system by evaluating integral between 
% sample points, and normalising by integral width (dt):
% [CJ(t)] * JLK = [C(t)] @ measurement times only
A(:,1) = (vCJ(second) - vCJ(first)) ./ dt;  
b = (vC(second) - vC(first)) ./ dt;
    
end