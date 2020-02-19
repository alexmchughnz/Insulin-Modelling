function [EGP] = find_EGP(EGP_init,y0,tol,SI,ksec,koffset)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%GC Constants


CNS = 0.3; %glucose consumption attributed to central nervous system (mmol·min-1)
VG = 12.2; % volume distribution of glucose (L)
nK = 0.06; %renal insulin clearance (min-1)
nL = 0.0324; %hepatic insulin clearance rate (min?1)
alphaI = 0.0017; %hepatic clearance saturation constant (L·mU-1)
xL = 0.67; %first pass constant as endogenous insulin is secreted into the portal vein 
nc = 0.032; %insulin degredation rate (min?1)
VI = 4; %volume distribution of insulin (L)
nI = 0.006; %27.56e-2/(VI); %the trans-endothelial diffusion rate between the plasma and interstitial compartments (min-1)
alphaG = 0.0154; %from lotz 2008, 0 for low dose injection, insulin binding saturation cinsulinonstant

EGP_1 = EGP_init;
EGP_new = 10;
maxiter = 100;
iter = 0;
while abs(EGP_new-EGP_1)>tol &&  iter<maxiter
    EGP_1 = EGP_new;
    dG_dt = @(t,y) - SI*y(1)*y(3) + (EGP_1 - CNS)/VG; %rate of change of plasma glucose concentration
    dI_dt = @(t,y) -nK*y(2) - nL*y(2)/(1+alphaI*y(2)) - nI*(y(2)-y(3)) + ((1-xL)*uen(y(1),ksec,koffset))/VI; %rate of change of plasma insulin concentration
    dQ_dt = @(t,y) nI*(y(2)-y(3)) - (nc * y(3))/(1+alphaG*y(3)); %rate of change of insulin in the interstitum

    odefun = @(t,y) [dG_dt(t,y); dI_dt(t,y);dQ_dt(t,y)];
    [GItime,GIglucose] = ode45(odefun, [0,1000], y0);
    
    vals = GIglucose(end,:);
    EGP_new = SI*vals(1)*vals(3)*VG +CNS;
    iter = iter + 1;
end
EGP = EGP_new;
iter
end

