%% function that performs a forward solve for C and Y, and calculates cpep secretion rate
function [output]=Cpep_Model(flag,Treal,Creal,S,val,timestep)
Y0=Creal(1)*val.ck1/val.ck2; %estimate for initial cpep concentration in peripheral compartment
switch(flag)
    case('ForwardSolve')
        [T_Cpep,C,Y_for] = Forward_Solve_Model(S,Creal(1),Y0,val,timestep);
        output.T_Cpep=T_Cpep; %time vector 
        output.C=C; %cpep concentration in central compartment
        output.Y_for=Y_for; %cpep concentration in peripheral compartment
    case('FindSecretion')
        [S,Y_Sec]=Find_Secretion_Model(Treal,Creal,Y0,val,timestep);
        output.S=S; %cpep secretion rate and time
        output.Y_Sec=Y_Sec; %cpep concentration in peripheral compartment
end
end

%% Forward Solve
function [T_Cpep,C,Y_for]=Forward_Solve_Model(S,C0,Y0,val,timestep)
format long
options=odeset('RelTol',1e-8,'AbsTol',1e-8);

ic=zeros(1,2); %Set up initial conditions
ic(1,1)=C0; %central cpep
ic(1,2)=Y0; %peripheral cpep

tfinal=[];
icfinal=[];
for i=1:length(timestep)-1 %loop through timesteps
    [t,ysol]=ode45(@Cpep_ODE,[timestep(i) timestep(i+1)],ic,options,val,S);
    ic=ysol(end,:)'; %sets up initial conditions for next timestep
    tfinal=[tfinal(1:end-1,:);t]; %put results onto the end of the last timestep
    icfinal=[icfinal(1:end-1,:);ysol];
end
T_Cpep=tfinal;
C=icfinal(:,1);
Y_for=icfinal(:,2);
end

%% ODEs for forward solve
function dCpepdt=Cpep_ODE(t,X,val,S)
Spp=interp1(S{1,1},S{1,2},t,'previous'); %creates interpolated fxn based on secretion rates and times

C=X(1,1);
Y=X(2,1);

dCpepdt=zeros(2,1);
dCpepdt(1)=(Spp)-((val.ck1+val.ck3)*C)+(val.ck2*Y);
dCpepdt(2)=(val.ck1*C)-(val.ck2*Y);
end

%% Secretion
function [S,Y_Sec]=Find_Secretion_Model(Treal,Creal,Y0,val,timestep)
format long

time=[timestep(1):1:timestep(end)]'; %create a time vector with 1 min intervals
Cpp=interp1(Treal,Creal,time); %use Treal and Creal to interpolate into 1 min vectors

Y_Sec=zeros(length(time),1); %Use C values and time points to calculate Y
Y_Sec(1,1)=Y0; %estimate for initial cpep concentration in peripheral compartment
for i=2:length(time)
    Y_Sec(i,1)=Y0.*exp(-val.ck2.*(time(i)+abs(time(1))))+trapz(time(1:i,1),val.ck1.*Cpp(1:i,1).*exp(-val.ck2.*(time(i,1)-time(1:i,1))));
end

Sec=zeros(length(time),1); %use Y values to calculate S
for i=1:length(time)-1
    Sec(i,1)=Cpp(i+1,1)-Cpp(i,1)+((val.ck1+val.ck3).*trapz(time(i:i+1,1),Cpp(i:i+1,1)))-(val.ck2.*(trapz(time(i:i+1,1),Y_Sec(i:i+1,1))));
end
Sec(length(time),1)=Sec(length(time)-1,1);

lb=1e-7;
Sec=max(Sec,lb);

S{1,1}=time; %time vector with 1 min intervals
S{1,2}=Sec; %S values for each minute interval
end