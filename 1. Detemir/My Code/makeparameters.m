% Adapted from "init_vars.m".

clear GI ID GC SC

%% Gastrointestinal (GI) Parameters
GI.k21 = 0.054;     %rate constant of grinding (min^-1)
GI.kAbs = 0.071;    %rate glucose absorbed into bloodstream (min^-1)
GI.kMax = 0.054;    %maximum value of kempt (min^-1)
GI.kMin = 0.006;    %minimum value of kempt (min^-1)
GI.f = 0.8;         %sclaing factor for incomplete absorption in first pass hepatic clearance

GI.DTot = 0.001;       %default amount of glucose in stomach - cannot be zero

%% Insulin Detemir (ID) Parameters
ID.k1 = 0.96;       %Detemir Unbind fraction
ID.k2 = 0.04;       %detemir bind fraction
ID.nK = 0.06;       %renal clearance
ID.nDL = 0.147;     %Detamir hepatic insulin clearance rate (min?1)
ID.nDI = 0.06;      %the detemir trans-endothelial diffusion rate between the plasma and interstitial compartments (min-1)
ID.nDC = 0.032;     %Detemir insulin degredation rate (min?1)
ID.C = 0.95;        %detemir bind/unbind rate

%% Glycaemic Control (GC) Parameters
GC.pg = 0.06;       % non insulin mediated uptake (min^-1)
GC.EGP = 0.96;      % endogenous glucose production (mmol/min)
GC.CNS = 0.3;       % glucose consumption attributed to central nervous system (mmol·min-1)
GC.VG = @(P) 0.18*P.mass;     % volume distribution of glucose (L)GC.nL = 0.15;       % hepatic insulin clearance rate (min?1)
GC.nK = 0.0644;     % NORMAL renal insulin clearance (min-1) - THIS GETS MODIFIED DURING THE FIRST DAY TO BE 1/2 OF THE VALUE, reverts to normal value after t_nk_change time has passed...
GC.alphaI = 0.0017; % hepatic clearance saturation constant (L·mU-1)
GC.xL = 0.67;       % first pass constant as endogenous insulin is secreted into the portal vein 
GC.nC = 0.032;      % insulin degredation rate (min?1)
GC.VI = @(P) 0.038*P.mass;   % volume distribution of insulin (L)
GC.nI = 0.006;      % 27.56e-2/(GC.VI); %the trans-endothelial diffusion rate between the plasma and interstitial compartments (min-1)
GC.alphaG = 0.0154; % from lotz 2008, 0 for low dose injection, insulin binding saturation cinsulinonstant

GC.uMin = 16.7;     % minimum endogenous insulin secretion (mU·min?1)
GC.uMax = 267;      % maximum endogenous insulin secretion (mU·min?1)
GC.nKChangeTime = 600;  % the time nk goes from nk/3 -> nk, also k3/3 -> k3
% GC.Uen = van_cauter(sys);   % solves Uen according to the van cauter model


%% Endogenous Insulin Secretion (SC) Parameters
SC.kdi = 0.0624;    % (min^-1)
SC.k2 = 0.0104;     %rate constant of diffusion from SC to interstitium(min^-1) 
SC.k3 = 0.0106;     %rate constant of insulin absorbed into plasma(min^-1)



save('parameters.mat', 'GI', 'ID', 'GC', 'SC')
disp('Parameters updated.')
clear