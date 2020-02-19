%Jacob Klenner - 99453904, jkl59@uclive.ac.nz
%Ben van Noorden - 15573361, bva22@uclive.ac.nz

%University of Canterbury Glycaemic Control ENME408 2018, M08

clear all, clc

%Index
%   GI = GastroIntestinal
%   SC = Subcutaneous
%   GC = Glycaemic Control
%   Data = Patient Data take

%Create Struct to Store Values
sys = struct('GI',[],'SC',[],'GC',[],'Data',[]);

%SBIT Trial Patient Data Analysis
%Patient Number
sys.Data.PtNo = 1; %Define Patient Number

%initialising ptient data
sys = PtDataRead(sys);

%initialising all variables
sys = init_vars(sys);
save(sprintf('patient%d.mat', sys.Data.PtNo))
%fitting insulin sensitivity values for the trial
sys = fit_SI(sys);

%GI Inital Conditions
q1 = 0.001; %initial amount of glucose in stomach compartment 1 (solid state)
q2 = 0; %initial amount of glucose in stomach compartment 2 (liquid state)
q3 = 0; %initial amount of glucose in gut compartment
Ra0 = 0;

%Insulin Detemir Initial Conditions
Isc = 0;
Dflocal = 0;
Dblocal = 0;
IDf = 0;
IDb = 0;
QDf = 0;
QDb = 0;

%GC Initial COnditions
G0 = sys.Data.bg3(find(sys.Data.bg3_time == sys.sim_start_t));
I0 = sys.Data.PlasmaI(1) * 1e-12 * 5808 / 33.7e-6 * 1000; %define initial plasma insulin concentration - convert from pmol/L to mU
Q0 = 9; %define initial interstitial insulin concentration

%Define Initial conditions
sys.IC.q0 = [q1; q2; q3; Isc; Dflocal; Dblocal; IDf; IDb; QDf; QDb; G0; I0; Q0]; %append initial conditions into vector form

%Solving the system
sys = eq_solve(sys);
sys.results.time = sys.results.time/(24*60) + sys.sim_start_t; %Convert model time into matlab datenum format to coincide with data

%Calculating RMS errors for the fits
t = 0 : sys.sim_max_t;
G_actual = ppval(sys.results.Gpp, t);
I_actual = ppval(sys.results.Ipp, t);
G_rms = sqrt(sum((G_actual - sys.results.G').^2)/sys.sim_max_t);
I_rms = sqrt(sum((I_actual - sys.results.I').^2)/sys.sim_max_t);

%Plotting
%Glucose fit
figure(1);
plot(sys.Data.bg3_time(2:end-1),sys.Data.bg3(2:end-1),'r*',sys.results.time,sys.results.G, 'k');   %,sys.Data.bg1_time,sys.Data.bg1,sys.Data.bg2_time,sys.Data.bg2)
legend('Blood Test','Model')
title('Plasma Glucose')
xlabel('Time')
ylabel('Plasma Glucose mmol/L')
datetick('x')
ylim([4 15])

%Insulin fit
figure(2);
plot(sys.Data.PlasmaI_time,sys.Data.PlasmaI,'r*',sys.results.time,sys.results.I*6.05+sys.results.IDf, 'k')
legend('Blood Test','Model')
title('Plasma Insulin')
xlabel('Time')
ylabel('Plamsa Insulin pmol/L')
datetick('x')
% ylim([0 2000])

%SI fit
figure(3);
plot(sys.results.time(1:2160), sys.GC.SI(1:2160), 'k')
title('Insulin Sensitivity')
xlabel('Time')
ylabel('S_I L/mU/min')
datetick('x')
% ylim([8e-4 12.5e-4])

figure(4);
plot(sys.results.time,sys.GC.Uen,'k')
title('Estimated Endogenous Insulin Secretion')
xlabel('Time')
ylabel('U_e_n mU/min')
datetick('x')
ylim([0 350])