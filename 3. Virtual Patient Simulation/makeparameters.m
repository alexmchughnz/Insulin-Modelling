clear C GI IN GC

%% Constants / Conversions (C)
C.MGlucose = 180.156;  % Glucose molar mass [g/mol]


%% Gastrointestinal (GI) Parameters
% Initial Values
GI.q0Sto1     = 35; % Stomach glucose in solid compartment [g->mmol]
GI.q0Sto2     = 0;  % Stomach glucose in liquid compartment [mmol]
GI.q0Gut      = 0;  % Gut glucose [mmol]

% Parameters
GI.k21 = 1;                % Rate constant of stomach grinding [1/min]
GI.D = 35/C.MGlucose*1000; % Amount of glucose ingested [g->mmol]
GI.kAbs = 0.205;           % Rate of intestinal glucose absorption [1/min]
GI.f = 0.8;                % Scaling factor for incomplete absorption [1]
GI.kMin = 0.013;           % Maximum value of kEmpt [1/min]
GI.kMax = 0.045;           % Minimum value of kEmpt [1/min]
GI.b = 0.85;               % Stomach emptying parameters [1]
GI.c = 0.25;               % ''

% Dependent Parameters
GI.Ra = @(qGut) GI.f * GI.kAbs * qGut;  % Rate of ingested glucose appearance [mmol/min]


%% Insulin Subcutaneous Injection (IN) Parameters
% Initial Values
IN.ISC0    = 2000;  % Subcutaneous insulin [mU]
IN.QLocal0 = 0;     % Insulin in local interstitium [mU]

% Parameters
IN.k2 = 0.0104;  % [1/min]
IN.k3 = 0.0613;  % [1/min]
IN.kdi = 0.0029;  % [1/min]
IN.IBolus = 2000;  % [mU]
IN.TBolus = 15;  % [min]


%% Glycaemic Control (GC)
% Initial Values
GC.G0 = 4.8;          % Blood glucose [mmol/L]
GC.I0 = 13.4185;      % Plasma insulin [mU/L] 
GC.Q0 = 7.4255;       % Interstitial insulin [mU/L]

% Parameters
GC.pg   = 0.06;       % Non insulin mediated uptake [1/min]
GC.GFast = 4.8;       % Fasting blood glucose level [mmol/L]
GC.alphaG = 0.0154;   % Insulin binding saturation parameter [L/mU]
GC.EGP = 0.612;       % Endogenous glucose production [mmol/min]
GC.CNS = 0.42;        % Central nervous system glucose uptake [mmol/min]
GC.VG = 13.3;         % Glucose volume of distribution [L]
GC.nK = 0.06;         % Renal insulin clearance [1/min]
GC.nL = 0.0324;       % Hepatic insulin clearance [1/min]
GC.alphaI = 0.0017;   % Hepatic clearance saturation constant [L/mU]
GC.nI = 0.2756;       % Trans-endothelial insulin diffusion [1/min]
GC.xL = 0.67;         % First-pass hepatic insulin clearance [1]
GC.VI = 4;            % Insulin volume of distribution [L]
GC.nC = 0.0324;       % Peripheral insulin degradation [1/min]
GC.uMin = 16.7;       % Minimum endogenous insulin secretion [mU/min]
GC.uMax = 267;        % Maximum endogenous insulin secretion [mU/min]
GC.kSec = 14.9;       % [mU*L/(mmol*min)]
GC.kOff = -49.9;      % [mU/min]

% Dependent Parameters
GC.uEx = @(QLocal) IN.k3 * QLocal;  % Exogenous insulin input [mmol/min]

%% Save
save('parameters.mat', 'C', 'GI', 'IN', 'GC')
disp('Parameters updated.')
clear