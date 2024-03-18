function P = EstimateInsulinSecretion(P)
% Estimates pancreatic insulin secretion rate (Uen) using a patient's 
% C-peptide data.
% Model from Van Cauter et al. (1992).
% Code based on DISST_nL_xL_test.m from J Ormsbee / P Docherty.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen

CP = P.parameters.CP;
GC = P.parameters.GC;
CONST = LoadConstants();

PrintStatusUpdate(P, "Estimating Uen..."); 

%% Setup
% Time of reading in sim [min]
% Concentraton of C-peptide [pmol/L]
[tCPep, vCPep] = GetData(P.data.CPep);

%% Interpolation
ppCPep = griddedInterpolant(tCPep, vCPep);

% Interpolate CPep values.
tArray = P.results.tArray;  % Time range [min]
dt = tArray(2) - tArray(1);
CPepArray = ppCPep(tArray);

%% Solving
% Initial conditions.
Y = zeros(size(tArray));
Y0 = CP.k1/CP.k2*CPepArray(1);
Y(1:2) = Y0;

% Calculate Y over time. Assumes dY/dt == 0(?)
for ii = 3:length(tArray)
    tSeg = tArray(1:ii-1);
    CPepSeg = CPepArray(1:ii-1);
    YSeg = Y(1:ii-1);
    Y(ii) = Y(1) + trapz(tSeg, CP.k1*CPepSeg - CP.k2*YSeg);
end

% Calculate endogenous secretion rate (Uen).
dvCPep = [diff(CPepArray)/dt; 0];
S = (dvCPep + (CP.k1+CP.k3).*CPepArray - CP.k2*Y);  % C-peptide secretion [(pmol/L)/min]
Uen = CONST.pmol2mU(S) * GC.VI;            % Endogenous insulin secretion [mU/min]

% Write value to patient struct.
P.results.Uen = Uen;

%% Plotting
MakePlots(P);

end

function MakePlots(P)
DP = DebugPlots().EstimateInsulinSecretion;

%% CPep
if DP.CPep
   MakeDebugPlot("C-peptide", P, DP);
   
   plt = plot(P.data.CPep.time, P.data.CPep.value, 'b*');
   plt.DisplayName = "Plasma Sample";
   
   xlabel("Time [min]")
   ylabel("Plasma C-Peptide, C [pmol/L]")
   
   legend
end

%% Uen
if DP.Uen
   MakeDebugPlot("Uen", P, DP);
   plot(P.results.tArray, P.results.Uen)
   
   xlabel("Time [min]")
   ylabel("Uen [mU/min]")
end

end
