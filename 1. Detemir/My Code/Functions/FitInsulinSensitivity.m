function P = FitInsulinSensitivity(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
%   GC  - glycaemic control parameter set
% OUTPUT:
%   P   - modified patient struct with SI

%% Setup


defaultSI = 10.8e-4;

% Define time range.
intervalTime = 360;  % Time per interval [min]
numIntervals = floor(P.trialLength/intervalTime);

% Find indices of trial times in CGM data.
iiStart = find(P.G{3}.time == P.trialTime(1));   % Trial start [index]
iiEnd   = find(P.G{3}.time == P.trialTime(end)); % Trial end [index]

% Get value/time arrays over trial period.
tG = minutes(P.G{3}.time(iiStart:iiEnd) - P.G{3}.time(iiStart));
vG = P.G{3}.value(iiStart:iiEnd);
tI = minutes(P.I.time - P.I.time(1));
vI = P.I.value;

ta = 1; 
tb = ta + intervalTime;

% Interpolate with a piecewise polynomial.
ppG = griddedInterpolant(tG, vG);
ppI = griddedInterpolant(tI, vI);

% Create SI arrays.
SI = ones(P.trialLength, 1) *defaultSI;


%% Compute

%STUB: Directly adapted from fit_SI.m

% ODE solver settings.
optionsShort = odeset('RelTol',1e-5,'AbsTol',1e-4,'MaxStep',0.2,'InitialStep',0.01);    % for first minute
optionsLong = odeset('RelTol',1e-5,'AbsTol',1e-4,'InitialStep',0.1);    % for subsequent minutes
ODEInitials = [0.001, 0, 0, 9, 0, 0];

for ii = 1 : numIntervals
    [TI1, Ints1] = ode45(@GIModelODE, (ta : ta+1)', ODEInitials, optionsShort, sys, Gpp, Ipp);
    ODEInitials = Ints1(end, :)';
    [TI2, Ints2] = ode45(@GIModelODE, (ta+1 : ta+intervalTime-1)', ODEInitials, optionsLong, sys, Gpp, Ipp);
    TI = [TI1(1, :); TI2];
    Ints = [Ints1(1, :); Ints2];


    A = Ints(:, 5);
    b = -diff(ppval(Gpp, [ta, tb]))' - Ints(:, 6);

    %update the conditions    
    ODEInitials = Ints(end, :);
    
    %Solve the linear system
    SI(ii) = A\b;
    
    ta = ta + intervalTime;
    temp_SI((ii-1)*intervalTime+1 : ii*intervalTime) = SI(ii);
    ta = ta + intervalTime;
    tb = tb + intervalTime;
end


P.SI(1 : length(temp_SI)) = temp_SI;

fprintf('SI fit successfully. SI ( 1e-4 ) = ')
disp(SI*1e4)
%\STUB

end


% Adapted from fit_SI/FAERIES_integrals.
function [dydt] = GIModelODE(t, y, ppG, ppI, GI, GC, P)

% Assign incoming variables.
qSto1  = y(1);
qSto2  = y(2);
qGut   = y(3);
Q      = y(4);

% Retrieve required derived values.
kEmpt = GetStomachEmptyingRate(t, (qSto1+qSto2), GI, P);
D = GetGlucoseDelivery(t, P);
G = ppG(t);  % Current blood glucose (interpolated) [mmol/L]
I = ppI(t);  % Current plasma insulin (interpolated) [mU/L]
Ra = GI.f * GI.kAbs * qGut;  % Rate of glucose appearance in plasma

% Gastrointestinal model DEs.
dqSto1 = D(t) - GI.k21*qSto1;
dqSto2 = GI.k21*qSto1;
dqGut  = kEmpt*qSto2 - GI.kAbs*qGut;

% Glycaemic control model DEs.
dQ     = GC.nI*(I - Q) - GC.nC * Q/(1+GC.alphaG*Q);
dGA    = G * Q/(1 + GC.alphaG*Q);
%dGb    = (Ra + EGP + infusion(sys, t) - sys.GC.CNS)/sys.GC.VG - sys.GC.pg*G;
dGb    =  - GC.pg*G + (Ra + EGP - GC.CNS)/GC.VG;
% NOTE: dG = -SI*dGA + dGb. SI is found with dGA\dGb (linear system).

% Pack up outgoing variables.
dydt = [dqSto1;
        dqSto2;
        dqGut;
        dQ;
        dGA;
        dGb];

end
