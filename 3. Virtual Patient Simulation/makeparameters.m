clear C GI IN SC

%% Constants / Conversions (C)
% All from J. L. Knopp's insulin conversion paper.
C.mIU2pmol = @(uIU) uIU * 6.0;    % Insulin [mIU]  -> Insulin [pmol]
C.pmol2mIU = @(pmol) pmol / 6.0;  % Insulin [pmol] -> Insulin [mIU]
C.MGlucose = 180.156;  % Glucose molar mass [g/mol]

C.IU18Factor = 18;  % TODO: Not sure what this yet! Something involving dL?

%% Gastrointestinal (GI) Parameters
GI.k21 = 1;     % Rate constant of stomach grinding [1/min]
GI.D = 25;   % Amount of glucose [mmol]
GI.kAbs = 0.205;    % Rate of intestinal glucose absorption [1/min]
GI.f = 0.8;         % Scaling factor for incomplete absorption [1]
GI.kMin = 0.013;    % Maximum value of kEmpt [1/min]
GI.kMax = 0.045;    % Minimum value of kEmpt [1/min]

GI.b = 0.85;        % Stomach emptying parameters [1]
GI.c = 0.25;        % ''
% GI.DTot0 = 0.001;   % Default amount of glucose in stomach [mmol]

%% Insulin (IN) Parameters
IN.k2 = 0.0104; % [1/min]
IN.k3 = 0.0613; % [1/min]
IN.kdi = 0.0029; % [1/min]
IN.IBolus = 2000; % [mU]
IN.T = 15; % [min]


%% Endogenous Insulin Secretion (SC) Parameters
F = 0.76;
CHalfLife1 = 4.95;       % Half-life of C-peptide in compartment 1 [min]
CHalflife2 = 32.4;         % Half-life of C-peptide in compartment 2 [min]
a = log(2)/CHalfLife1;
b = log(2)/CHalflife2;
% NOTE: Currently unused in EstimateInsulinSecretion!
SC.k2 = F*(b-a) + a;     % Rate constants
SC.k3 = a*b/(2*SC.k2); % NOTE: Original code had no factor of 1/2; PDD's thesis does.
SC.k1 = a + b - SC.k2 - SC.k3;

%% Glycaemic Control (GC) Parameters
GC.pg   = 0.06;                       % Non insulin mediated uptake [1/min]
GC.GFast = 4.8; %[mmol/L]
GC.SI = 10.8e-4; %[L/(mU*min)]
GC.alphaG = 0.0154; %[L/mU]
GC.EGP = 0.612;                        % Endogenous glucose production [mmol/min]
GC.CNS = 0.42;                         % Central nervous system glucose uptake [mmol/min]
GC.VG = 13.3;             % Glucose volume of distribution [L]
GC.nK = 0.06;                       % Renal insulin clearance [1/min]
GC.nL = 0.0324;                         % Hepatic insulin clearance [1/min]
GC.alphaI = 0.0017;                   % Hepatic clearance saturation constant [L/mU]
GC.nI = 0.2756;                        % Trans-endothelial insulin diffusion [1/min]
GC.xL = 0.67;                         % First-pass hepatic insulin clearance [1]
GC.VI = 4;            % Insulin volume of distribution [L]
GC.nC = 0.0324;                        % Peripheral insulin degradation [1/min]
GC.uMin = 16.7;                       % Minimum endogenous insulin secretion [mU/min]
GC.uMax = 267;                        % Maximum endogenous insulin secretion [mU/min]
GC.kSec = 14.9;     % [mU*L/(mmol*min)]
GC.kOff = -49.9;     % [mU/min]


%%
save('parameters.mat', 'C', 'GI', 'IN', 'GC')
disp('Parameters updated.')
clear