global GC

% Parameters
GC.pg   = 0.004;                       % Non insulin mediated uptake [1/min]
GC.EGP = 0.96;                        % Endogenous glucose production [mmol/min]
GC.CNS = 0.3;                         % Central nervous system glucose uptake [mmol/min]
GC.nK = 0.0644;                       % Renal insulin clearance [1/min]
GC.alphaI = 0.0017;                   % Hepatic clearance saturation constant [L/mU]
GC.alphaG = 0.0154;                   % Saturation of insulin binding to cells [L/mU]

GC.uMin = 16.7;                       % Minimum endogenous insulin secretion [mU/min]
GC.uMax = 267;                        % Maximum endogenous insulin secretion [mU/min]


% Dependent Parameters
global SC
GC.VI = 4;                    % Insulin volume of distribution [L]
GC.VQ = (SC.k1/SC.k2)*GC.VI;  % Interstitial volume of distribution [L]
GC.VG = 1.2*(GC.VI + GC.VQ);  % Glucose volume of distribution [L]

GC.nI = GC.VQ * SC.k2;                 % Trans-endothelial insulin diffusion [1/min]
gamma = 0.5;
GC.nC = (GC.nI/GC.VQ) * (1/gamma - 1); % Peripheral insulin degradation [1/min]