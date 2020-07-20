function P = EstimateInsulinSecretion(P)
% Estimates pancreatic insulin secretion rate (Uen) using a patient's 
% C-peptide data.
% Model from Van Cauter et al. (1992).
% Code based on van_cauter.m.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen


global C SC GC
global DEBUGPLOTS

%% Setup
% Time of reading in sim [min]
% Concentraton of C-peptide [pmol/L]
[tCPep, vCPep] = GetSimTime(P, P.data.CPep);
vCPep = vCPep * GC.VI; % Total C-peptides [pmol]

% Rate constants.
k1 = SC.k1;   
k2 = SC.k2;   
k3 = SC.k3;

%% Interpolation
ppCPep = griddedInterpolant(tCPep, vCPep);

% Make 1-minute spaced time vector, and interpolate CPep values.
t = P.results.tArray;  % Time range [min]
dt = t(2) - t(1);
vCPep = ppCPep(t);

%% Solving
% Calculate Y over time. Assumes dY/dt == 0;
Y = zeros(size(t));
Y(1) = k1(1)/k2*vCPep(1);
Y = exp(-k2*t) .* (Y(1) + cumtrapz(t, exp(k2*t).*k1.*vCPep));

% Calculate endogenous secretion rate (Uen).
S = [diff(vCPep); 0]/dt + (k1 + k3).*vCPep - k2*Y;  % C-peptide secretion [pmol/min]
Uen = C.pmol2mU(S);                                 % Endogenous insulin secretion [mU/min]

% Write value to patient struct.
P.results.Uen = Uen;
fprintf('P%d: Estimated Uen.\n', P.patientNum); 

%% Debug Plots
DP = DEBUGPLOTS.EstimateInsulinSecretion;
if DP.Uen
   MakeDebugPlot(P, DP);
   plot(t, Uen)
   title(sprintf("%s: Uen", P.patientCode))
   
   xlabel("Time [min]")
   ylabel("Uen [mU/min]")
end

end

