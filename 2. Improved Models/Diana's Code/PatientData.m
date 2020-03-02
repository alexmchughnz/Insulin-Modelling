%list of parameters used in GC_Model_Solver
%from Holder-Pearson et al 2018
dt_insulin = 0.1;
dt_dex = 0.5; %30 sec
%Lui_Patient_1b
lui.patient="Patient 1b";
lui.delta_T_ins=dt_insulin;
lui.delta_T_dex=dt_dex;
lui.t_fp=[-22;-1;10;20;26;33;41;50;61;70;79;91;100;121];
lui.G_fp=[5.2;5;6.2;8;10.6;10.5;10.3;9.3;6.9;6.7;6;5.5;4.9;4.8];

% lui.G_iv=[4.5;5.9;8;11.1;9.4;7.1;5.5;4.5];
% lui.C_iv=[681;973;1370;1950;1810;1570;1250;925];
% lui.I_iv=[26;107;187;228;160;101;47;34]/6;
% lui.t_iv=[-32;11;20;41;50;61;91;121];

lui.G_iv=[4.5;5.9;8;11.1;7.1;5.5;4.5];
lui.I_iv=[26;107;187;228;101;47;34]/6;
lui.t_iv=[-32;11;20;41;61;91;121];

lui.C_iv=[681;1625;973;1370;1950;1570;1250;925]; %"ideal" cpep
% lui.C_iv=[681;1950;973;1370;1950;1570;1250;925]; %max cpep
lui.tc_iv=[-32;5;11;20;41;61;91;121];

lui.D_bolus_size=(35*1000)/(180);
lui.D_bolus_time=0;
lui.Dinput{1,1}=[lui.t_iv(1);lui.D_bolus_time;lui.D_bolus_time+lui.delta_T_dex];
lui.Dinput{1,2}=[0;lui.D_bolus_size/lui.delta_T_dex;0];
lui.I_bolus_size_1=1000;
lui.I_bolus_size_2=1800;
lui.I_bolus_time_1=15;
lui.I_bolus_time_2=17;
lui.Ibolusinput{1,1}=[lui.t_iv(1);lui.I_bolus_time_1;lui.I_bolus_time_1+lui.delta_T_ins;lui.I_bolus_time_2;lui.I_bolus_time_2+lui.delta_T_ins];
lui.Ibolusinput{1,2}=[0;lui.I_bolus_size_1/lui.delta_T_ins;0;lui.I_bolus_size_2/lui.delta_T_ins;0];
lui.timestep=[lui.Ibolusinput{1,1};lui.Dinput{1,1};lui.t_iv(end)];
lui.timestep = sort(lui.timestep);
lui.timestep = unique(lui.timestep);

save('Patient_01b','lui')
clear lui

%Lui_Patient_1c
lui.patient="Patient 1c";
lui.delta_T_ins=dt_insulin;
lui.delta_T_dex=dt_dex;
lui.t_fp=[-20;0;10;16;20;26;30;41;51;62;71;81;93;100;123];
lui.G_fp=[5.2;5.2;6.5;7.8;8.5;9.9;9.9;9.5;7.8;6.5;5.9;5.2;4.9;4.7;4.2];

% lui.G_iv=[5.2;5.2;8.5;9.5;7.8;6.5;4.9;4.2];
% lui.C_iv=[718;717;1470;1940;1840;1590;856;783];
% lui.I_iv=[42;41;193;236;165;94;10;24]/6;
% lui.t_iv=[-20;0;20;41;51;62;93;123];

lui.G_iv=[5.2;5.2;8.5;9.5;6.5;4.9;4.2];
lui.I_iv=[42;41;193;236;94;10;24]/6;
lui.t_iv=[-20;0;20;41;62;93;123];

lui.C_iv=[718;717;1530;1470;1940;1590;856;783]; %"ideal" cpep
% lui.C_iv=[718;717;1940;1470;1940;1590;856;783]; %max cpep
lui.tc_iv=[-20;0;5;20;41;62;93;123];

lui.D_bolus_size=(35*1000)/(180);
lui.D_bolus_time=0;
lui.Dinput{1,1}=[lui.t_iv(1);lui.D_bolus_time;lui.D_bolus_time+lui.delta_T_dex];
lui.Dinput{1,2}=[0;lui.D_bolus_size/lui.delta_T_dex;0];
lui.I_bolus_size=1800;
lui.I_bolus_time=15;
lui.Ibolusinput{1,1}=[lui.t_iv(1);lui.I_bolus_time;lui.I_bolus_time+lui.delta_T_ins];
lui.Ibolusinput{1,2}=[0;lui.I_bolus_size/lui.delta_T_ins;0];
lui.timestep=[lui.Ibolusinput{1,1};lui.Dinput{1,1};lui.t_iv(end)];
lui.timestep = sort(lui.timestep);
lui.timestep = unique(lui.timestep);
save('Patient_01c','lui')
clear lui

%Lui_Patient_2a
lui.patient="Patient 2a";
lui.delta_T_ins=dt_insulin;
lui.delta_T_dex=dt_dex;
lui.t_fp=[-47;-8;12;23;27;32;41;50;60;70;80;91;100;119];
lui.G_fp=[5.1;5.2;4.7;6.2;6.4;7.2;8;7.5;8.5;8.8;8.3;7.2;7.4;6.8];

% lui.G_iv=[5;4.9;6.5;8.4;8.2;8.7;7.5;6.2];
% lui.C_iv=[544;389;793;1330;1430;1780;1720;1370];
% lui.I_iv=[62;27;126;203;212;270;188;113]/6;
% lui.t_iv=[-47;-8;23;41;50;60;91;120];

lui.G_iv=[5;4.9;6.5;8.4;8.7;7.5;6.2];
lui.I_iv=[62;27;126;203;270;188;113]/6;
lui.t_iv=[-47;-8;23;41;60;91;120];

lui.C_iv=[544;389;890;793;1330;1780;1720;1370]; %"ideal" cpep
% lui.C_iv=[544;389;1780;793;1330;1780;1720;1370]; %max cpep
lui.tc_iv=[-47;-8;5;23;41;60;91;120];

lui.D_bolus_size=(35*1000)/(180);
lui.D_bolus_time=0;
lui.Dinput{1,1}=[lui.t_iv(1);lui.D_bolus_time;lui.D_bolus_time+lui.delta_T_dex];
lui.Dinput{1,2}=[0;lui.D_bolus_size/lui.delta_T_dex;0];
lui.I_bolus_size=1000;
lui.I_bolus_time=15;
lui.Ibolusinput{1,1}=[lui.t_iv(1);lui.I_bolus_time;lui.I_bolus_time+lui.delta_T_ins];
lui.Ibolusinput{1,2}=[0;lui.I_bolus_size/lui.delta_T_ins;0];
lui.timestep=[lui.Ibolusinput{1,1};lui.Dinput{1,1};lui.t_iv(end)];
lui.timestep = sort(lui.timestep);
lui.timestep = unique(lui.timestep);
save('Patient_02a','lui')
clear lui

%Lui_Patient_4a
lui.patient="Patient 4a";
lui.delta_T_ins=dt_insulin;
lui.delta_T_dex=dt_dex;
lui.t_fp=[-6;1;13;17;25;31;39;50;61;70;81;96;99;120];
lui.G_fp=[4.7;4.9;4.7;5.5;4.9;5.2;5.7;6.8;6.4;7.4;7.3;6.7;6.6;5];

% lui.G_iv=[4.6;4.7;4.6;4.7;5.9;7;6.5;7;5.2];
% lui.C_iv=[344;347;350;346;590;748;784;1100;735];
% lui.I_iv=[23;21;23;36;125;144;125;148;64]/6;
% lui.t_iv=[0;1;13;23;39;50;61;91;120];

lui.G_iv=[4.6;4.7;4.6;4.7;5.9;6.5;7;5.2];
lui.I_iv=[23;21;23;36;125;125;148;64]/6;
lui.t_iv=[0;1;13;23;39;61;91;120];

lui.C_iv=[344;347;300;350;346;590;784;1100;735]; %"ideal" cpep
% lui.C_iv=[344;347;1100;350;346;590;784;1100;735]; %max cpep
lui.tc_iv=[0;1;5;13;23;39;61;91;120];

lui.D_bolus_size=(35*1000)/(180);
lui.D_bolus_time=lui.t_iv(2);
lui.Dinput{1,1}=[lui.t_iv(1);lui.D_bolus_time;lui.D_bolus_time+lui.delta_T_dex];
lui.Dinput{1,2}=[0;lui.D_bolus_size/lui.delta_T_dex;0];
lui.I_bolus_size=2000;
lui.I_bolus_time=15;
lui.Ibolusinput{1,1}=[lui.t_iv(1);lui.I_bolus_time;lui.I_bolus_time+lui.delta_T_ins];
lui.Ibolusinput{1,2}=[0;lui.I_bolus_size/lui.delta_T_ins;0];
lui.timestep=[lui.Ibolusinput{1,1};lui.Dinput{1,1};lui.t_iv(end)];
lui.timestep = sort(lui.timestep);
lui.timestep = unique(lui.timestep);
save('Patient_04a','lui')
clear lui

%Lui_Patient_05a
lui.patient="Patient 5a";
lui.delta_T_ins=dt_insulin;
lui.delta_T_dex=dt_dex;
lui.t_fp=[-2;8;15;21;28;43;55;67;71;81;96;102;126];
lui.G_fp=[5.9;6.5;7;6.2;8.5;9.5;10.1;9;8.7;7.9;5.1;5.9;4.2];

% lui.G_iv=[6.2;6.5;10;10.2;9.2;5.1;4];
% lui.C_iv=[608;1040;2330;2810;2840;1340;680];
% lui.I_iv=[48;157;524;631;538;123;42]/6;
% lui.t_iv=[-2;18;48;51;62;73;122];

lui.G_iv=[6.2;6.5;10;9.2;5.1;4];
lui.I_iv=[48;157;524;538;123;42]/6;
lui.t_iv=[-2;18;48;62;73;122];

lui.C_iv=[608;1800;1040;2330;2840;1340;680]; %"ideal" cpep
% lui.C_iv=[608;2840;1040;2330;2840;1340;680]; %max cpep
lui.tc_iv=[-2;5;18;48;62;73;122];

lui.D_bolus_size=(35*1000)/(180);
lui.D_bolus_time=0;
lui.Dinput{1,1}=[lui.t_iv(1);lui.D_bolus_time;lui.D_bolus_time+lui.delta_T_dex];
lui.Dinput{1,2}=[0;lui.D_bolus_size/lui.delta_T_dex;0];
lui.I_bolus_size=1800;
lui.I_bolus_time=15;
lui.Ibolusinput{1,1}=[lui.t_iv(1);lui.I_bolus_time;lui.I_bolus_time+lui.delta_T_ins];
lui.Ibolusinput{1,2}=[0;lui.I_bolus_size/lui.delta_T_ins;0];
lui.timestep=[lui.Ibolusinput{1,1};lui.Dinput{1,1};lui.t_iv(end)];
lui.timestep = sort(lui.timestep);
lui.timestep = unique(lui.timestep);
save('Patient_05a','lui')
clear lui

%Lui_Patient_16a
lui.patient="Patient 16a";
lui.delta_T_ins=dt_insulin;
lui.delta_T_dex=dt_dex;
lui.t_fp=[-16;0;11;15;20;26;30;41;50;60;71;81;90;100;116];
lui.G_fp=[5.2;5.3;6.4;6.6;6;7.2;7;7.6;6.3;5.3;5.9;5.5;5.2;5.2;4.5];

% lui.G_iv=[5.15;5.11;6.43;7.1;6.67;5.54;4.99;4.67];
% lui.C_iv=[769;1230;2040;2080;1880;1430;1010;725];
% lui.I_iv=[44;110;217;137;105;54;24;30]/6;
% lui.t_iv=[-16;0;20;39;50;60;90;118];

lui.G_iv=[5.15;5.11;6.43;7.1;5.54;4.99;4.67];
lui.I_iv=[44;110;217;137;54;24;30]/6;
lui.t_iv=[-16;0;20;39;60;90;118];

lui.C_iv=[769;1230;2425;2040;2080;1430;1010;725]; %"ideal" cpep
% lui.C_iv=[769;1230;2080;2040;2080;1430;1010;725]; %max cpep
lui.tc_iv=[-16;0;5;20;39;60;90;118];

lui.D_bolus_size=(35*1000)/(180);
lui.D_bolus_time=0;
lui.Dinput{1,1}=[lui.t_iv(1);lui.D_bolus_time;lui.D_bolus_time+lui.delta_T_dex];
lui.Dinput{1,2}=[0;lui.D_bolus_size/lui.delta_T_dex;0];
lui.I_bolus_size=1800;
lui.I_bolus_time=15;
lui.Ibolusinput{1,1}=[lui.t_iv(1);lui.I_bolus_time;lui.I_bolus_time+lui.delta_T_ins];
lui.Ibolusinput{1,2}=[0;lui.I_bolus_size/lui.delta_T_ins;0];
lui.timestep=[lui.Ibolusinput{1,1};lui.Dinput{1,1};lui.t_iv(end)];
lui.timestep = sort(lui.timestep);
lui.timestep = unique(lui.timestep);
save('Patient_16a','lui')
clear lui
