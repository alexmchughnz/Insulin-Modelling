function [dY] = GCModelODEFast(t, Y, P, Y0)
% ODE for GC model. Use with ode45.
% Requires P1(t) and P2(t) - must be run AFTER GI model.
% Requires QDF(t) - must be run AFTER ID model.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [G; I; Q] at time == t
%   P   - patient struct
%   Y0  - initial conditions of states
% OUTPUT:
%   dY  - derivatives of states at time == t

global GC

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
GInfusion = P.data.GInfusion(n);    % Glucose infusion rate [mmol/min]
GFast     = P.data.GFast(t);        % Fasting glucose [mmol/L]
GBolus    = P.data.GBolus(t);       % Glucose bolus [mmol/min]
IBolus    = P.data.IBolus(t);       % Insulin bolus [mU/min]

% Patient dependent.
d2 = P.results.d2;
VG = GC.VG(P);
VQ = GC.VQ(P);
nI = GC.nI(P);
nC = GC.nC(P);
nK = GC.nK(P);

%% Computation
dG  = -GC.pg*(G-GFast) ...
          - SI*(G*Q - GFast*Q0)/(1 + GC.alphaG*Q) ...
          + d2/VG*P2 + (GInfusion + GBolus)/VG;      
dI  = -nK*I - nL*I/(1 + GC.alphaI*I) - nI/GC.VI*(I-Q) ...
          + Uen*(1 - xL)/GC.VI + IBolus/GC.VI;
dQ  = nI/VQ*(I-Q) - nC*Q;

%% Output
dY = [dG;
      dI;
      dQ];

end
