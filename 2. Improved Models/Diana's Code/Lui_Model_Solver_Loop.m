clear,clc
for i=1:6
    
    val=initiate_val_struct(); %struct with parameters used for calculations
    
    if i==1
       load('Patient_01b.mat')
        d2=-log(0.5)./[2:2:20];
    elseif i==2
        load('Patient_01c.mat')
        d2=-log(0.5)./[2:2:20];
    elseif i==3
        load('Patient_02a.mat')
        d2=-log(0.5)./[30:2:70];
    elseif i==4
        load('Patient_04a.mat')
        d2=-log(0.5)./[70:2:110];
    elseif i==5
        load('Patient_05a.mat')
        d2=-log(0.5)./[2:2:20];
    elseif i==6
        load('Patient_16a.mat')
        d2=-log(0.5)./[2:2:20];
    end
    
    output=Cpep_Model('FindSecretion',lui.tc_iv,lui.C_iv,[],val,lui.timestep);
    val.S=output.S;
    val.Y_Sec=output.Y_Sec;
    val.Uen=val.S{1,2}*val.Vp/6;

    output=Cpep_Model('ForwardSolve',[],lui.C_iv,val.S,val,lui.timestep);
    val.T_Cpep=output.T_Cpep;
    val.C=output.C;
    val.Y_for=output.Y_for;

    [nL,xL]=nLxL_Model(lui.t_iv,lui.I_iv,val,lui.timestep,lui.Ibolusinput);
    val.nL=nL;
    val.xL=xL;

    [SI,d2]=Grid_Search_Model(d2,val,lui.t_iv,lui.G_iv,lui.I_iv,lui.timestep,lui.Dinput,lui.Ibolusinput);
    val.SI=SI;
    val.d2=d2;

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

    str=lui.patient+": SI = "+num2str(val.SI)+", d2 = "+num2str(val.d2)+", nL = "+num2str(val.nL)+", xL = "+num2str(val.xL)+", ErrG = "+num2str(val.errg)+", Erri = "+num2str(val.erri);

    figure()
    grid on
%     subplot(4,1,1)
    subplot(5,1,1)
    plot(lui.t_iv,lui.G_iv,'ro-',lui.t_fp,lui.G_fp,'rx',val.T_GIQ,val.G,'b-')
    title(str)
    ylabel('BGL (mmol/L)')
    xlim([-50 140])
    ylim([4 12])
    legend('IV Measured','FP Measured','Model')

%     subplot(4,1,2)
    subplot(5,1,2)
    plot(lui.t_iv,lui.I_iv,'ro-',val.T_GIQ,val.I,'b-',val.T_GIQ,val.Q,'g-')
    ylabel('Insulin (mU/L)')
    xlim([-50 140])
    ylim([0 100])
    legend('Measured I','Model I','Model Q')

%     subplot(4,1,3)
    subplot(5,1,3)
    plot(lui.tc_iv,lui.C_iv,'ro-',val.T_Cpep,val.C,'b-')
    ylabel('Central C-Pep Concentration (pmol/L)')
    xlim([-50 140])
    ylim([0 2500])
    legend('Measured','Model')

%     subplot(4,1,4)
    subplot(5,1,4)
    plot(val.S{1,1},val.Uen,'b-')
    ylabel('Uen (mU/min)')
    xlim([-50 140])
    ylim([0 250])
    legend('Model')
    
    subplot (5,1,5)
    plot(lui.t_iv,lui.I_iv-I)
    ylabel('Measured I - Simulated I')
    xlim([-50 140])
    
    xlabel('Time (min)')
    
    pause(0.1)
end