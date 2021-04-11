function [A, b, IFunc, QFunc] = AssembleIIntegralSystem(P, I, Q)

CONST = LoadConstants();
GC = P.parameters.GC;



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
% Consider dI/dt = k + cx*(1-xL) - kI*I - cn*nL - kIQ*(I-Q):
cx = Uen/GC.VI;
cn = I./(1 + GC.alphaI*I);

k = IInput/GC.VI;
kI = GC.nK;
kIQ = GC.nI./GC.VI;

% Also consider dQ/dt = -cQ*Q + cI*I:
cQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q - easier to use
cI = GC.nI/GC.VQ;  % Constant term coefficent of I - easier to use

%% Integrate I Equation
% I(t) - I(t0) = int{k} + int{cx}*(1-xL) - kI*int{I} - int{cn}*nL - kIQ*int{I-Q}
% Defining CN = -int{cn} and CX = int{cx}
% CN*nL + CX*(1-xL) = I(t) - I(t0) - int{k} + kI*int{I} + kIQ*int{I-Q} := C
CN = -cumtrapz(tArray, cn);
CX = cumtrapz(tArray, cx);

intk = cumtrapz(tArray, k);
intI = kI*cumtrapz(tArray, I);
intIQ = kIQ*cumtrapz(tArray, I-Q);

I0 = I(1) * ones(size(I));
RHS = [I -I0 -intk intI intIQ];
C = sum(RHS, CONST.ROWWISE);

%% Make Minute-Wise Q and I Functions
% I(t) = I(t0) + int{k} + CX*(1-xL) - kI*int{I} + CN*nL - kIQ*int{I-Q}
IFunc = @(nL, xL, I, Q) I(1) + cumtrapz(tArray, k) + CX*(1-xL) - kI*cumtrapz(tArray, I) ...
    + CN*nL - kIQ*cumtrapz(tArray, I-Q);

% Q(t) = Q(t0) - cQ*int{Q} + cI*int{I}
QFunc = @(I, Q) Q(1) - cQ*cumtrapz(tArray, Q) + cI*cumtrapz(tArray, I);

%% Assemble MLR System
% Extract values at measurement points.
vCN = CN(iiMeas);
vCX = CX(iiMeas);
vC = C(iiMeas);

% Find time width of each integral wedge. 
first = 1:numel(vC) - 1;
second = 2:numel(vC);

t1 = tMeas(first);
t2 = tMeas(second);
dt = t2 - t1; 

% Assemble MLR system by evaluating integral between 
% sample points, and normalising by integral width (dt):
% [CN(t) CX(t)] * (nL; 1-xL) = [C(t)] @ measurement times only
A(:,1) = (vCN(second) - vCN(first)) ./ dt;  
A(:,2) = (vCX(second) - vCX(first)) ./ dt;
b = (vC(second) - vC(first)) ./ dt;
    
end