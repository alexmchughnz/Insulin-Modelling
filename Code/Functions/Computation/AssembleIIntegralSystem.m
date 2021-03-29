function [A, b, LHS, RHS] = AssembleIIntegralSystem(P, tMeas, I, Q)

ROWWISE = 2;

GC = P.parameters.GC;

%% Setup
% tArray = P.results.tArray;
% t0 = tArray(1);

iiMeas = GetTimeIndex(tMeas, P.results.tArray);

% Plasma Insulin
vI = I(iiMeas);

% Interstitial Insulin
vQ = Q(iiMeas);

% Endogenous Secretion
vUen = P.results.Uen(iiMeas);

% Exogenous Insulin
vIInput = GetPlasmaInsulinInput(tMeas, P);

%% Get Coefficients
% Consider dI/dt = k + cx*(1-xL) - kI*I - cn*nL - kIQ*(I-Q):
cx = vUen/GC.VI;
cn = vI./(1 + GC.alphaI*vI);

k = vIInput/GC.VI;
kI = GC.nK;
kIQ = GC.nI./GC.VI;

% Also consider dQ/dt = -cQ*Q + cI*I:
cQ = GC.nC + GC.nI/GC.VQ; % Constant term coefficent of Q - easier to use
cI = GC.nI/GC.VQ;  % Constant term coefficent of I - easier to use


%% Integrate I Equation
% I(t) - I(t0) = int{k} + int{cx}*(1-xL) - kI*int{I} - int{cn}*nL - kIQ*int{I-Q}
% Defining CN = -int{cn} and CX = int{cx}
% CN*nL + CX*(1-xL) = I(t) - I(t0) - int{k} + kI*int{I} + kIQ*int{I-Q} := C
CN = -cumtrapz(tMeas, cn);
CX = cumtrapz(tMeas, cx);
intk = cumtrapz(tMeas, k);
intI = kI*cumtrapz(tMeas, vI);
intIQ = kIQ*cumtrapz(tMeas, vI-vQ);

vI0 = I(1) * ones(size(vI));

LHS = [CN CX];
RHS = [vI -vI0 -intk intI intIQ];
C = sum(RHS, ROWWISE);

%% Assemble MLR System
first = 1:numel(vI) - 1;
second = 2:numel(vI);

t1 = tMeas(first);
t2 = tMeas(second);
dt = t2 - t1;  % Time width of each wedge.

% Area of each integral wedge for each component.   
vCN = CN(second) - CN(first);
vCX = CX(second) - CX(first);
vC  = C(second) - C(first);

% Assembling MLR system, integrating between sample points, and
% normalising by integral width (dt):
% [CN(t) CX(t)] * (nL; 1-xL) = [C(t)]
A(:,1) = vCN ./ dt;  
A(:,2) = vCX ./ dt;
b = vC ./ dt;

end