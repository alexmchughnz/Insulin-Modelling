function [sys] = init_vars(sys)
%Contains ODE's and solver for a subcut insulin, gastrointestinal glucose
%and glycaeimic control model
% 
%Define Model Constants
%GI Constants
sys.GI.kabs = 0.071;    %rate glucose absorbed into bloodstream (min^-1)
sys.GI.f = 0.8;         %sclaing factor for incomplete absorption in first pass hepatic clearance
sys.GI.kmax = 0.054;    %maximum value of kempt (min^-1)
sys.GI.kmin = 0.006;    %minimum value of kempt (min^-1)
sys.GI.k21 = 0.054;     %rate constant of grinding (min^-1)
sys.GI.D = 0.001;       %default amount of glucose in stomach - cannot be zero

%SC Constants
sys.SC.kdi = 0.0624;    % (min^-1)
sys.SC.k2 = 0.0104;     %rate constant of diffusion from SC to interstitium(min^-1) 
sys.SC.k3 = 0.0106;     %rate constant of insulin absorbed into plasma(min^-1)

%Detemir Constants
sys.ID.C = 0.95;        %detemir bind/unbind rate
sys.ID.k1 = 0.96;       %Detemir Unbind fraction
sys.ID.k2 = 0.04;       %detemir bind fraction
sys.ID.nDI = 0.06;      %the detemir trans-endothelial diffusion rate between the plasma and interstitial compartments (min-1)
sys.ID.nDc = 0.032;     %Detemir insulin degredation rate (min?1)
sys.ID.nK = 0.06;       %renal clearance
sys.ID.nDL = 0.147;     %Detamir hepatic insulin clearance rate (min?1)

%GC Constants
sys.GC.pg = 0.06;       % non insulin mediated uptake (min^-1)
sys.GC.CNS = 0.3;       % glucose consumption attributed to central nervous system (mmol·min-1)
sys.GC.VG = 0.18*sys.Data.pt_mass;     % volume distribution of glucose (L)
sys.GC.nK = 0.0644;     % NORMAL renal insulin clearance (min-1) - THIS GETS MODIFIED DURING THE FIRST DAY TO BE 1/2 OF THE VALUE, reverts to normal value after t_nk_change time has passed...
sys.GC.nL = 0.15;       % hepatic insulin clearance rate (min?1)
sys.GC.alphaI = 0.0017; % hepatic clearance saturation constant (L·mU-1)
sys.GC.xL = 0.67;       % first pass constant as endogenous insulin is secreted into the portal vein 
sys.GC.nc = 0.032;      % insulin degredation rate (min?1)
sys.GC.VI = 0.038*sys.Data.pt_mass;   % volume distribution of insulin (L)
sys.GC.nI = 0.006;      % 27.56e-2/(sys.GC.VI); %the trans-endothelial diffusion rate between the plasma and interstitial compartments (min-1)
sys.GC.umin = 16.7;     % minimum endogenous insulin secretion (mU·min?1)
sys.GC.umax = 267;      % maximum endogenous insulin secretion (mU·min?1)
sys.GC.alphaG = 0.0154; % from lotz 2008, 0 for low dose injection, insulin binding saturation cinsulinonstant
sys.GC.t_nk_change = 600;  % the time nk goes from nk/3 -> nk, also k3/3 -> k3
sys.GC.Uen = van_cauter(sys);   % solves Uen according to the van cauter model
sys.GC.EGP = 0.96;      % endogenous glucose production (mmol/min)

%EGP SS Solve - bit of a black box, don't really know how or why it does
%this but it seems to work?
% sys.GC.EGP = 1; %EGP initial guess
% sys.GC.EGP = suggestEGP(sys,0.0001,100); %find_EGP(EGP_init,y_infin,tol,SI,ksec,koffset);

%set EGP to constant initially - work on this when you get all the other
%stuff working

disp('Initialised system parameters')
end