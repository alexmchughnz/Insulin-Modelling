function [dY] = GCModelODE(t, Y, P, Y0, ppG, ppI)
% ODE for GC model. Use with ode45.
% Note that G terms are split into those multiplying SI (GA) and those
% added to it (Gb). This allows fitting SI from forward sim.
% Requires P1(t) and P2(t) - must be run AFTER GI model.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [G; I; Q; GA; Gb] at time == t
%   P   - patient struct
%   Y0  - initial conditions of states
%   ppG - (optional) profile of glucose over time to use
%   ppI - (optional) profile of insulin over time to use
% OUTPUT:
%   dY  - derivatives of states at time == t

GC = P.parameters.GC;

%% Input
G  = Y(1);  % [mmol/L]
I  = Y(2);  % [mU]
Q  = Y(3);  % [mU]
Q0 = Y0(3); % [mU]

%% Variables
% Time dependent.
n      = GetTimeIndex(t, P.results.tArray);  % Index of current timestep.
Uen    = P.results.Uen(n);       % [mU/min]
P2     = P.results.P2(n);        % [mmol]
GFast  = P.data.GFast(t);        % Fasting glucose [mmol/L]

% Patient dependent.
SI = P.results.SI;        % [L/mU/min]
d2 = P.results.d2;
nL = P.results.nL;        % [1/min]
xL = P.results.xL;        % [1]

% Trial dependent.
if P.data.GDelivery == "intravenous"
    GInput = P.data.GBolus(t) + P.data.GInfusion(n);  % [mmol/min]
else
    GInput = P.data.GInfusion(n);  % [mmol/min]
end


Qtot = Q;
Qtot0 = Q0;
if P.data.IType == "detemir"
    % Include interstitial Detemir.
    Qtot = Qtot + P.results.QDF(n);   % [mU/L]
    Qtot0 = Qtot0 + P.results.QDF(1);  % [mU/L]
end

IInput = GetPlasmaInsulinInput(t, P);

%% Computation
dGA = -(G*Qtot - GFast*Qtot0)/(1 + GC.alphaG*Qtot);   % Insulin-mediated uptake.

dGb =GInput/GC.VG...          % Exogenous IV glucose input.
    + d2/GC.VG * P2 ...       % Endogenous input from gut.
    + GC.EGP/GC.VG ...        % Endogenous production.
    - GC.pg * (G-GFast) ...   % Non insulin-mediated uptake.
    - GC.CNS/GC.VG;           % Central nervous system uptake.

dG = dGb + SI*dGA;

dI  = IInput/GC.VI ...              % Exogenous input (IV or subcut).
    + Uen * (1 - xL)/GC.VI ...      % Endogenous input.
    - GC.nK * I ...                 % Renal clearance.
    - nL * I/(1 + GC.alphaI*I) ...  % Hepatic clearance.
    - GC.nI/GC.VI * (I-Q);          % Transfer to Q compartment.


dQ  = GC.nI/GC.VQ * (I-Q) ...  % Transfer from I compartment.
    - GC.nC * Qtot;            % Peripheral degradation.



%% Output
dY = [dG;
    dI;
    dQ;
    dGA;
    dGb];

end
