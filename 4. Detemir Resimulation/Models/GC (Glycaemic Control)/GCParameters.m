% Parameters
GC.pg   = 0.06;                       % Non insulin mediated uptake [1/min]
GC.EGP = 0.96;                        % Endogenous glucose production [mmol/min]
GC.CNS = 0.3;                         % Central nervous system glucose uptake [mmol/min]
GC.nK = 0.0644;                       % Renal insulin clearance [1/min]
GC.nL = 0.15;                         % Hepatic insulin clearance [1/min]
GC.alphaI = 0.0017;                   % Hepatic clearance saturation constant [L/mU]
GC.xL = 0.67;                         % First-pass hepatic insulin clearance [1]
GC.nC = 0.032;                        % Peripheral insulin degradation [1/min]
GC.nI = 0.006;                        % Trans-endothelial insulin diffusion [1/min]
GC.alphaG = 0.0154;                   % Saturation of insulin binding to cells [L/mU]

GC.uMin = 16.7;                       % Minimum endogenous insulin secretion [mU/min]
GC.uMax = 267;                        % Maximum endogenous insulin secretion [mU/min]

% Dependent Parameters
GC.VG = @(P) 0.18*P.mass;             % Glucose volume of distribution [L]
GC.VI = @(P) 0.038*P.mass;            % Insulin volume of distribution [L]
GC.VQ = @(P, SC) (SC.k1/SC.k2)*GC.VI(P);  % Interstitial volume of distribution [L]