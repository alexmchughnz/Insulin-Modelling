global GC

% Parameters
GC.pg   = 0.004;                       % Non insulin mediated uptake [1/min]
GC.EGP = 0.96;                        % Endogenous glucose production [mmol/min]
GC.CNS = 0.3;                         % Central nervous system glucose uptake [mmol/min]
GC.alphaI = 0.0017;                   % Hepatic clearance saturation constant [L/mU]
GC.alphaG = 0.0154;                   % Saturation of insulin binding to cells [L/mU]


% Dependent Parameters
GC.VI = 4;                    % Insulin volume of distribution [L]
GC.VQ = @(P) (P.data.k1/P.data.k2)*GC.VI;  % Interstitial volume of distribution [L]
GC.VG = @(P) 1.2*(GC.VI + GC.VQ(P));  % Glucose volume of distribution [L]

GC.nI = @(P) GC.VQ(P) * P.data.k2;                 % Trans-endothelial insulin diffusion [1/min]
gamma = 0.5;
GC.nC = @(P) (GC.nI(P)/GC.VQ(P)) * (1/gamma - 1); % Peripheral insulin degradation [1/min]

GC.nK = @(P) P.data.k3;                       % Renal insulin clearance [1/min]
