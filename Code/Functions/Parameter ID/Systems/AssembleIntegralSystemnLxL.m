function [A, b, IFunc, QFunc, paramNames] = AssembleIntegralSystemnLxL(P, I, Q)

CONST = Constants();
GC = P.parameters.GC;

paramNames = ["nL", "xL"];

%% Setup
tArray = P.results.tArray;
tMeas = P.data.I.time;
iiMeas = SearchArray(tMeas, tArray);

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
Uex = P.results.Uex(P);

%% Get Coefficients
% Consider dI/dt = k + cx*1 - cx*xL - cn*nL - kI*I - kIQ*(I-Q):
cx = Uen/GC.VI;
cn = I./(1 + GC.alphaI*I);

k = Uex/GC.VI;
kI = GC.nK;
kIQ = GC.nI./GC.VI;

% Also consider dQ/dt = -cQ*Q + cI*I:
cQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q - easier to use
cI = GC.nI/GC.VQ;  % Constant term coefficent of I - easier to use

%% Integrate I Equation
% I(t) - I(t0) = int{k} + int{cx}*1 - int{cx}*xL - int{cn}*nL - kI*int{I} - kIQ*int{I-Q}
% Defining CN = int{cn} and CX = int{cx}
% CN*nL + CX*xL = -I(t) + I(t0) + int{k} + CX - kI*int{I} - kIQ*int{I-Q} := C
CN = cumtrapz(tArray, cn);
CX = cumtrapz(tArray, cx);

intkTerm = cumtrapz(tArray, k);
intITerm = kI*cumtrapz(tArray, I);
intIQTerm = kIQ*cumtrapz(tArray, I-Q);

I0 = I(1) * ones(size(I));
RHS = [-I +I0 +intkTerm +CX -intITerm -intIQTerm ];
C = sum(RHS, CONST.COLUMNDIM);

%% Make Minute-Wise Q and I Functions
% I(t) = I(t0) + int{k} + CX - kI*int{I} - kIQ*int{I-Q} - CN*nL - CX*xL
IFunc = @(nL, xL, I, Q) I(1) + intkTerm + CX - intITerm - intIQTerm - CN*nL - CX*xL;

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
% [CN(t) CX(t)] * (nL; xL) = [C(t)] @ measurement times only
A(:,1) = (vCN(second) - vCN(first)) ./ dt;  
A(:,2) = (vCX(second) - vCX(first)) ./ dt;
b = (vC(second) - vC(first)) ./ dt;
    
end