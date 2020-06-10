function P = PDEstimateInsulinSecretion(P)
% Estimates pancreatic insulin secretion rate (Uen) using a patient's 
% C-peptide data.
% Model from Van Cauter et al. (1992).
% Code based on van_cauter.m.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen


global C GC SC

%% Setup
% Time and data arrays.
t = minutes(P.G{3}.time(1:end-1) - P.G{3}.time(1));  % Time or reading [min]
CPep0 = P.CPep.value(1);

UB = SC.k3 * CPep0 * GC.VI(P); %Basal insulin secretion [pmol/min]
G = P.G{3}.value(1:end-1);
dG = diff(P.G{3}.value);
GB = P.GFast(t);

%% Solving    
gainP = 1;
gainD = 3 * (dG > 0);

Uen = UB + gainP*(G - GB) + gainD.*dG;

%% Interpolation
ppUen = griddedInterpolant(t, Uen);

% Make 1-minute spaced time vector, and interpolate Uen values.
tInterpolated = (0 : P.simDuration-1)';  % Time range,  [min]
UenInterpolated = ppUen(tInterpolated);

% Write value to patient struct.
P.Uen.value = UenInterpolated;
P.Uen.time  = tInterpolated;
fprintf('P%d: Estimated Uen.\n', P.patientNum); 

end

