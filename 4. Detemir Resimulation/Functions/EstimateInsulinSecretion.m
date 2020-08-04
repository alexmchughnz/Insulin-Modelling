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
% vCPep = vCPep * GC.VI; % Total C-peptides [pmol]

% Rate constants.
k1 = SC.k1;   
k2 = SC.k2;   
k3 = SC.k3;

%% Interpolation
ppCPep = griddedInterpolant(tCPep, vCPep);

% Make 1-minute spaced time vector, and interpolate CPep values.
tArray = [P.data.simTime(1) : P.data.simTime(end)]';  % Time range [min]
dt = tArray(2) - tArray(1);
vCPep = ppCPep(tArray);

%% Solving
% Initial conditions.
Y = zeros(size(tArray));
Y0 = k1/k2*vCPep(1);
Y(1:2) = Y0;

% Calculate Y over time. Assumes dY/dt == 0(?)
for ii = 3:length(tArray)
    tSeg = tArray(1:ii-1);
    CPepSeg = vCPep(1:ii-1);
    YSeg = Y(1:ii-1);
    Y(ii) = Y(1) + trapz(tSeg, k1*CPepSeg - k2*YSeg);
end

% Calculate endogenous secretion rate (Uen).
dvCPep = [diff(vCPep)/dt; 0];
S = (dvCPep + (k1+k3).*vCPep - k2*Y);  % C-peptide secretion [(pmol/L)/min]
Uen = C.pmol2mU(S) * GC.VI;            % Endogenous insulin secretion [mU/min]

% Write value to patient struct.
P.results.Uen = Uen;
fprintf('P%d: Estimated Uen.\n', P.patientNum); 

%% Debug Plots
DP = DEBUGPLOTS.EstimateInsulinSecretion;
if DP.Uen
   MakeDebugPlot(P, DP);
   plot(tArray, Uen)
   title(sprintf("%s: Uen", P.patientCode))
   
   xlabel("Time [min]")
   ylabel("Uen [mU/min]")
end

end

