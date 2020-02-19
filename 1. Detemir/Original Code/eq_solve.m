function [sys] = eq_solve(sys)
%the function that solves the system of ode's



% sys.SC.Ibolus=0;
%Define Equations
%Gastrointestinal model
dq1_dt = @(t,y) - sys.GI.k21*y(1) + eating(sys, t);% + carbohydrates(sys, t);% + glucose(sys, t); %rate of change of glucose in stomach solid state (note: Dirac term replaced with initial condition)
dq2_dt = @(t,y) - kempt(t,y(1),y(2),sys)*y(2) + sys.GI.k21*y(1); %rate of change of glucose in liquid solid state
dq3_dt = @(t,y) - sys.GI.kabs*y(3) + kempt(t,y(1),y(2),sys)*y(2); %rate of change of glucose in gut

%Insulin Detemir model
dIsc_dt = @(t,y) -sys.SC.k2*y(4) + utotal(t,sys); %rate of change of subcut insulin
dDflocal_dt = @(t,y) -sys.SC.k3*y(5) + sys.SC.k2*y(4) - sys.SC.kdi*y(5) - sys.ID.C*(sys.ID.k1*y(5) - sys.ID.k2*y(6)); %rate of change of insulin in the local interstitium
dDblocal_dt = @(t,y) sys.ID.C*(sys.ID.k1*y(5) - sys.ID.k2*y(6));
dIDf_dt = @(t,y) -sys.ID.nDL*y(7) - sys.ID.nDI*(y(7)-y(9)) - sys.ID.C*(sys.ID.k1*y(7)-sys.ID.k2*y(8)) + sys.SC.k3*y(5)/sys.GC.VI - sys.ID.nK*y(7);
dIDb_dt = @(t,y) sys.ID.C*(sys.ID.k1*y(7) - sys.ID.k2*y(8));
dQDf_dt = @(t,y) sys.ID.nDI*(y(7)-y(9)) - sys.ID.nDc*y(9) - sys.ID.C*(sys.ID.k1*y(9) - sys.ID.k2*y(10));
dQDb_dt = @(t,y) sys.ID.C*(sys.ID.k1*y(9) - sys.ID.k2*y(10));

%Glycaemic control model
dG_dt = @(t,y) - sys.GC.pg*(y(11) - fasting_BG(sys,t)) - get_SI(sys,t)*y(11)*(y(13)+y(9))/(1 + sys.GC.alphaG*(y(13)+y(9))) + (sys.GI.f*sys.GI.kabs*y(3) + sys.GC.EGP + infusion(sys, t) - sys.GC.CNS)/sys.GC.VG; %rate of change of plasma glucose concentration
dI_dt = @(t,y) - nk(sys, t)*y(12) - sys.GC.nL*y(12)/(1 + sys.GC.alphaI*y(12)) - sys.GC.nI*(y(12) - y(13)) + ((1 - sys.GC.xL)*uen(sys, t))/sys.GC.VI; %rate of change of plasma insulin concentration
dQ_dt = @(t,y) sys.GC.nI*(y(12) - y(13)) - (sys.GC.nc*y(13))/(1 + sys.GC.alphaG*y(13)); %rate of change of insulin in the interstitium

%Append equations to vector form
odefun = @(t,y) [dq1_dt(t,y); dq2_dt(t,y); dq3_dt(t,y); dIsc_dt(t,y); dDflocal_dt(t,y); dDblocal_dt(t,y); dIDf_dt(t,y); dIDb_dt(t,y); dQDf_dt(t,y); dQDb_dt(t,y); dG_dt(t,y); dI_dt(t,y); dQ_dt(t,y)];

% %Solve 1 (from time 0 to time T when SC insulin is injected)
% tspan1 = 0 : 0.5 : sys.SC.T;
% [GItime1,GIglucose1] = ode45(odefun, tspan1, sys.IC.q0);
% 
% %Solve 2 (from time T to time 120)
% tspan2 = sys.SC.T : 0.5 : sys.max_t;
% sys.IC.q02 = GIglucose1(end,:); %take initial conditions from end of first solve
% sys.IC.q02(4) = sys.SC.Ibolus; %substitute in Ibolus (replaces dirac term)
% 
% [GItime2,GIglucose2] = ode45(odefun, tspan2, sys.IC.q02);

%Solve 1 (from time 0 to time T when SC insulin is injected)
tspan1 = (0 : sys.sim_max_t);
[GItime,GIglucose] = ode45(odefun, tspan1, sys.IC.q0);


%append results
sys.results.time = GItime;
sys.results.q1 = GIglucose(:,1);
sys.results.q2 = GIglucose(:,2);
sys.results.q3 = GIglucose(:,3);
sys.results.Isc = GIglucose(:,4);
sys.results.Dflocal = GIglucose(:,5);
sys.results.Dblocal = GIglucose(:,6);
sys.results.IDf = GIglucose(:,7)*18;
sys.results.IDb = GIglucose(:,8)*18;
sys.results.QDf = GIglucose(:,9)*18;
sys.results.QDb = GIglucose(:,10)*18;
sys.results.G = GIglucose(:,11);
sys.results.I = GIglucose(:,12);
sys.results.Q = GIglucose(:,13);

sys.results.serumDet = sys.results.IDf+sys.results.IDb;

% sys.results.Ra = GIglucose(:,11);
% sys.results.Ins_plasma = GIglucose(:,12);
% sys.results.Ins_lost = GIglucose(:,13);
end