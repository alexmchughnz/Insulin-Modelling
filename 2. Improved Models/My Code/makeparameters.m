% Adapted from "init_vars.m".

clear C GI ID GC SC

%% Constants / Conversions (C)
% From J. L. Knopp's insulin conversion paper.
C.mIU2pmol = @(uIU) uIU * 6.0;    % Insulin [mIU]  -> Insulin [pmol]
C.pmol2mIU = @(pmol) pmol / 6.0;  % Insulin [pmol] -> Insulin [mIU]
C.MGlucose = 180.156;  % Glucose molar mass [g/mol]

C.IU18Factor = 18;  % TODO: Not sure what this yet! Something involving dL?

%% Gastrointestinal (GI) Parameters
P1HalfLife = 20;    % Half-life of glucose in stomach (P1) [min]
P2HalfLife = 60;    % Half-life of glucose in gut (P2) [min]
GI.d1 = log(2)/P1HalfLife;
GI.d2 = log(2)/P2HalfLife;  % NOTE: dummy value, needs parameter ID.

GI.k21 = 0.054;     % Rate constant of stomach grinding [1/min]
GI.kAbs = 0.071;    % rate glucose absorbed into bloodstream [1/min]
GI.b = 0.69;        % Stomach emptying parameters [1]
GI.c = 0.17;        % ''
GI.kMax = 0.054;    % Maximum value of kEmpt [1/min]
GI.kMin = 0.006;    % Minimum value of kEmpt [1/min]
GI.f = 0.8;         % Scaling factor for incomplete absorption [1]

GI.DTot0 = 0.001;   % Default amount of glucose in stomach [mmol]

%% Insulin Detemir (ID) Parameters
ID.ka   = 0.0078;   % Hexameric Detemir dissociation rate [1/min]
ID.kd1  = 0.96;     % Detemir unbind fraction [1/min]
ID.kd2  = 0.04;     % Detemir bind fraction [1/min]
ID.kdi  = 0.0594;   % Detemir interstitial degradation [1/min]
ID.kb   = 0.0181;   % Unbound I diffusion from local to blood [1/min]
% NOTE: nK doesn't match Klenner paper (nK = 0).
ID.nK   = 0.06;     % Detemir renal clearance [1/min]
ID.nDL  = 0.147;    % Detamir hepatic insulin clearance rate [1/min]
ID.nDI  = 0.06;     % Detemir trans-endothelial diffusion rate [1/min]
ID.nDC  = 0.032;    % Detemir insulin degredation rate [1/min]

%% Endogenous Insulin Secretion (SC) Parameters
F = 0.76;
CHalfLife1 = 4.95;       % Half-life of C-peptide in compartment 1 [min]
CHalflife2 = 32.4;         % Half-life of C-peptide in compartment 2 [min]
a = log(2)/CHalfLife1;
b = log(2)/CHalflife2;
SC.k2 = F*(b-a) + a;     % Rate constants
SC.k3 = a*b/(2*SC.k2); % NOTE: Original code had no factor of 1/2; PDD's thesis does.
SC.k1 = a + b - SC.k2 - SC.k3;

%% Glycaemic Control (GC) Parameters
GC.pg   = 0.06;                       % Non insulin mediated uptake [1/min]
GC.EGP = 0.96;                        % Endogenous glucose production [mmol/min]
GC.CNS = 0.3;                         % Central nervous system glucose uptake [mmol/min]
GC.VG = @(P) 0.18*P.mass;             % Glucose volume of distribution [L]
GC.nK = 0.0644;                       % Renal insulin clearance [1/min]
GC.nL = 0.15;                         % Hepatic insulin clearance [1/min]
GC.alphaI = 0.0017;                   % Hepatic clearance saturation constant [L/mU]
GC.xL = 0.67;                         % First-pass hepatic insulin clearance [1]
GC.nC = 0.032;                        % Peripheral insulin degradation [1/min]
GC.VI = @(P) 0.038*P.mass;            % Insulin volume of distribution [L]
GC.nI = 0.006;                        % Trans-endothelial insulin diffusion [1/min]
GC.alphaG = 0.0154;                   % Saturation of insulin binding to cells [L/mU]
GC.VQ = @(P) (SC.k1/SC.k2)*GC.VI(P);  % Interstitial volume of distribution [L]

GC.uMin = 16.7;                       % Minimum endogenous insulin secretion [mU/min]
GC.uMax = 267;                        % Maximum endogenous insulin secretion [mU/min]

%%
save('parameters.mat', 'C', 'GI', 'ID', 'SC', 'GC')
disp('Parameters updated.')
clear