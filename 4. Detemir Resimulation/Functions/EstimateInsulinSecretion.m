function P = EstimateInsulinSecretion(P)
% Estimates pancreatic insulin secretion rate (Uen) using a patient's 
% C-peptide data.
% Model from Van Cauter et al. (1992).
% Code based on DISST_nL_xL_test.m from J Ormsbee / P Docherty.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen


global C GC
global DEBUGPLOTS

%% Setup
% Time of reading in sim [min]
% Concentraton of C-peptide [pmol/L]
[tCPep, vCPep] = GetSimTime(P, P.data.CPep);

% Rate constants.
k1 = P.data.k1;   
k2 = P.data.k2;   
k3 = P.data.k3;

%% Interpolation
ppCPep = griddedInterpolant(tCPep, vCPep);

% Make 1-minute spaced time vector, and interpolate CPep values.
tArray = [P.data.simTime(1) : P.data.simTime(end)-1]';  % Minute-wise time range [min]
dt = tArray(2) - tArray(1);
CPepArray = ppCPep(tArray);

%% Solving
% Initial conditions.
Y = zeros(size(tArray));
Y0 = k1/k2*CPepArray(1);
Y(1:2) = Y0;

% Calculate Y over time. Assumes dY/dt == 0(?)
for ii = 3:length(tArray)
    tSeg = tArray(1:ii-1);
    CPepSeg = CPepArray(1:ii-1);
    YSeg = Y(1:ii-1);
    Y(ii) = Y(1) + trapz(tSeg, k1*CPepSeg - k2*YSeg);
end

% Calculate endogenous secretion rate (Uen).
dvCPep = [diff(CPepArray)/dt; 0];
S = (dvCPep + (k1+k3).*CPepArray - k2*Y);  % C-peptide secretion [(pmol/L)/min]
Uen = C.pmol2mU(S) * GC.VI;            % Endogenous insulin secretion [mU/min]

% Expand to original time, and write value to patient struct.
dt = P.results.tArray(2) - P.results.tArray(1);
P.results.Uen = repelem(Uen, 1/dt, 1);
fprintf('P%d: Estimated Uen.\n', P.patientNum); 

%% Debug Plots
DP = DEBUGPLOTS.EstimateInsulinSecretion;
if DP.CPep
   MakeDebugPlot(P, DP);
   hold on
   
   plt = plot(tCPep, vCPep, 'r*');
   plt.DisplayName = "Data";
   plt = plot(tArray, CPepArray, 'b-');
   plt.DisplayName = "Interpolation";
   
   title(sprintf("%s: C-peptide samples", P.patientCode))
   
   xlabel("Time [min]")
   ylabel("C-peptide [pmol/L]")
   
   legend
end
if DP.Uen
   MakeDebugPlot(P, DP);
   plot(tArray, Uen)
   title(sprintf("%s: Uen", P.patientCode))
   
   xlabel("Time [min]")
   ylabel("Uen [mU/min]")
end

end

