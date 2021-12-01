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

% Apply smoothing.
k = 2;
centralWeight = 0.6;

filter = ones(2*k+1, 1) * (1-centralWeight)/(2*k);
filter(k+1) = centralWeight;


smoothUen = conv2(Uen, filter);
smoothUen = smoothUen(k+1:end-k);

% Write value to patient struct.
P.results.Uen = smoothUen;

%% Plotting
plotvars.Uen = Uen;
plotvars.smoothUen = smoothUen;

MakePlots(P, plotvars);

end

function MakePlots(P, plotvars)
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
   plt = plot(P.results.tArray, plotvars.Uen, '.');
   plt.DisplayName = "Raw Uen";

   plt = plot(P.results.tArray, plotvars.smoothUen, '.');
   plt.DisplayName = "Smooth Uen";

   xlabel("Time [min]")
   ylabel("Uen [mU/min]")
end

end
