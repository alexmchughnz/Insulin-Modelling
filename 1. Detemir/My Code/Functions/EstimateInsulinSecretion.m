function P = EstimateInsulinSecretion(P)
% Estimates pancreatic insulin secretion rate (Uen) using a patient's 
% C-peptide data.
% Model from Van Cauter et al. (1992).
% Code based on van_cauter.m.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen


global C GC
%% Data and Parameters
% Hard coded rate constants.
% NOTE: why not the SC model parameter values?
k1 = 0.0478;    %+/- 0.019
k2 = 0.0516;    %+/- 0.013

% Volume of insulin distribution.
% NOTE: why not the SC model parameter values?
V = 5; %+/- 0.83

% Time and data arrays.
t = minutes(P.CPep.time - P.CPep.time(1));  % Time or reading [min]
CPep = P.CPep.value * V;  % Amount of C-peptide [pmol]


%% Interpolation
ppCPep = griddedInterpolant(t, CPep);

% Make 1-minute spaced time vector, and interpolate CPep values.
dt = 1;
t = (0 : dt : t(end))';  % Time range,  [min]
CPep = ppCPep(t);


%% Solving

% Set up time dependant parameters.
k1 = k1 * ones(size(t));
%k3 = GC.nK(1:n)';
k3 = GC.nK(1) * ones(size(t));  % NOTE: original code has this constant...?

% Calculate Y over time. Assumes dY/dt == 0;
Y = zeros(size(t));
Y(1) = k1(1)/k2*CPep(1);
Y = exp(-k2*t).*(Y(1) + cumtrapz(t,exp(k2*t).*k1.*CPep));

% Calculate endogenous secretion rate (Uen).
Uen = [diff(CPep); 0]/dt + (k1 + k3).*CPep - k2*Y;  % [pmol/min]
Uen = C.pmol2mIU(Uen);                              % [mU/min]

% Write value to patient struct.
P.Uen.value = Uen;
P.Uen.time  = t;
fprintf('P%d: Estimated Uen.\n', P.patientNum); 

end

