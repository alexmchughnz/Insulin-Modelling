function [dY] = GCModelODE(t, Y, P, Y0)
% ODE for GC model. Use with ode45.
% Note that G terms are split into those multiplying SI (GA) and those
% added to it (Gb). This allows fitting SI from forward sim.
% Requires P1(t) and P2(t) - must be run AFTER GI model.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [G; I; Q; GA; Gb] at time == t
%   P   - patient struct
%   Y0  - initial conditions of states
% OUTPUT:
%   dY  - derivatives of states at time == t


GC = P.parameters.GC;
SC = P.parameters.SC;

%% Input
G  = Y(1);  % [mmol/L]
I  = Y(2);  % [mU]
Q  = Y(3);  % [mU]
Q0 = Y0(3); % [mU]

%% Variables
% Time dependent.
n = GetTimeIndex(t, P.results.tArray);  % Index of current timestep.
SI        = P.results.SI(n);        % [L/mU/min]
Uen       = P.results.Uen(n);       % [mU/min]
P2        = P.results.P2(n);        % [mmol]
nL        = P.results.nL(n);        % [1/min]
xL        = P.results.xL(n);        % [1]
GFast     = P.data.GFast(t);        % Fasting glucose [mmol/L]

% Patient dependent.
d2 = P.results.d2;

% Trial dependent.
if P.data.GDelivery == "intravenous"
    GInput = P.data.GBolus(t) + P.data.GInfusion(n);  % [mmol/min]
else
    GInput = P.data.GInfusion(n);  % [mmol/min]
end

if P.data.IType == "detemir"
    % Include interstitial Detemir.
    Q0 = Q0 + P.results.QDF(1);  % [mU/L]
    Q = Q + P.results.QDF(n);   % [mU/L]
end

IInput = GetPlasmaInsulinInput(t, P);

%% Computation
dGA = -(G*Q - GFast*Q0)/(1 + GC.alphaG*Q);   % Insulin-mediated uptake.

dGb = - GC.pg*(G-GFast) ...                  % Non-insulin-mediated uptake.
    + d2*P2/GC.VG ...                           % Endogenous input from gut.
    + GInput/GC.VG;                             % Exogenous IV glucose input.

dG = SI*dGA + dGb;

dI  = - GC.nK*I ...                  % Renal clearance.
    - nL*I/(1 + GC.alphaI*I) ...  % Liver clearance.
    - GC.nI/GC.VI*(I-Q) ...          % Transfer to Q compartment.
    + Uen*(1 - xL)/GC.VI ...      % Endogenous input.
    + IInput/GC.VI;               % Exogenous input (IV or subcut).

dQ  = GC.nI/GC.VQ*(I-Q) ...  % Transfer from I compartment.
    - GC.nC*Q;            % Peripheral degradation.



%% Output
dY = [dG;
    dI;
    dQ;
    dGA;
    dGb];

end
