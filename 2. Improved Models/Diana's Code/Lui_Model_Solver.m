clear,clc

d2=-log(0.5)./[10:10:90];

val=initiate_val_struct(); %struct with parameters used for calculations
load('Patient_01b.mat') %load in patient data from Lui's trials

% Calculates cpep secretion rate based on measured cpep values
output=Cpep_Model('FindSecretion',lui.tc_iv,lui.C_iv,[],val,lui.timestep);
val.S=output.S;
val.Y_Sec=output.Y_Sec;
val.Uen=val.S{1,2}*val.Vp/6;

% Calculates central and peripheral cpep concentrations
output=Cpep_Model('ForwardSolve',[],lui.C_iv,val.S,val,lui.timestep);
val.T_Cpep=output.T_Cpep;
val.C=output.C;
val.Y_for=output.Y_for;

% Calculates nL and xL
[nL,xL]=nLxL_Model(lui.t_iv,lui.I_iv,val,lui.timestep,lui.Ibolusinput);
val.nL=nL;
val.xL=xL;

[SI,d2]=Grid_Search_Model(d2,val,lui.t_iv,lui.G_iv,lui.I_iv,lui.timestep,lui.Dinput,lui.Ibolusinput);
val.SI=SI;
val.d2=d2;

% Calculates G, I, and Q
[output,val.T_GIQ]=GIQ_Model(val,lui.G_iv,lui.I_iv,lui.timestep,lui.Dinput,lui.Ibolusinput);
val.P1=output(:,1);
val.P2=output(:,2);
val.Isc=output(:,3);
val.Qlocal=output(:,4);
val.G=output(:,5);
val.I=output(:,6);
val.Q=output(:,7);

G=interp1(val.T_GIQ,val.G,lui.t_iv);
I=interp1(val.T_GIQ,val.I,lui.t_iv);
errg=((abs(lui.G_iv-G))./lui.G_iv)*100;
erri=((abs(lui.I_iv-I))./lui.I_iv)*100;
val.errg=median(errg);
val.erri=median(erri);

str=lui.patient+": SI = "+num2str(val.SI)+", da = "+num2str(val.da)+", d2 = "+num2str(val.d2)+", nL = "+num2str(val.nL)+", xL = "+num2str(val.xL)+", ErrG = "+num2str(val.errg)+", Erri = "+num2str(val.erri);

figure()
subplot(4,1,1)
plot(lui.t_iv,lui.G_iv,'ro-',lui.t_fp,lui.G_fp,'rx',val.T_GIQ,val.G,'b-')
title(str)
ylabel('BGL (mmol/L)')
xlim([-50 140])
ylim([4 12])
legend('IV Measured','FP Measured','Model')

subplot(4,1,2)
plot(lui.t_iv,lui.I_iv,'ro-',val.T_GIQ,val.I,'b-',val.T_GIQ,val.Q,'g-')
ylabel('Insulin (mU/L)')
xlim([-50 140])
ylim([0 100])
legend('Measured I','Model I','Model Q')

subplot(4,1,3)
plot(lui.tc_iv,lui.C_iv,'ro-',val.T_Cpep,val.C,'b-')
ylabel('Central C-Pep Concentration (pmol/L)')
xlim([-50 140])
ylim([0 2500])
legend('Measured','Model')

subplot(4,1,4)
plot(val.S{1,1},val.Uen,'b-')
ylabel('Uen (mU/min)')
xlim([-50 140])
ylim([0 250])
legend('Model')
xlabel('Time (min)')
