%% function to solve for SI given BGL measurements
function [SI]=SI_Model(val,Treal,Greal,Ireal,timestep,Dinput,Iinput)
format long

Gpp=interp1(Treal,Greal,'linear','pp'); %fit a piecewise fxn to the BGL measurements

ic=zeros(1,8); %initial conditions of ODEs
ic(1,1)=0; %P1
ic(1,2)=0; %P2
ic(1,3)=0; %Isc
ic(1,4)=0; %Qlocal
ic(1,5)=Ireal(1); %I
ic(1,6)=Ireal(1)/2; %Q
ic(1,7)=0; %A (used to solve for SI)
ic(1,8)=0; %b (used to solve for SI)

options=odeset('RelTol',1e-8,'AbsTol',1e-8);
tfinal=[];
icfinal=[]; 
for i=1:length(timestep)-1 %loop through timesteps
    [t,ysol]=ode45(@SI_Model_ODE,[timestep(i) timestep(i+1)],ic,options,val,Greal,Ireal,Gpp,Dinput,Iinput);
    ic=ysol(end,:)'; %sets up initial conditions for next timestep
    tfinal=[tfinal(1:end-1,:);t]; %put results onto the end of the last timestep
    icfinal=[icfinal(1:end-1,:);ysol];
end

A=zeros(1,1); %Construct linear system A*sI=b
b=zeros(1,1);
A(1,1)=icfinal(find(tfinal<=timestep(end),1,'last'),7)-icfinal(find(tfinal<=timestep(1),1,'last'),7);
b(1,1)=-(ppval(Gpp,timestep(end))-ppval(Gpp,timestep(1)))+(icfinal(find(tfinal<=timestep(end),1,'last'),8)-icfinal(find(tfinal<=timestep(1),1,'last'),8));
SI=A\b;
lb=1e-7;
SI=max(SI,lb);
end

%% function that sets up differential equations to solve for SI given BGL measurements
function dSIdt=SI_Model_ODE(t,X,val,Greal,Ireal,Gpp,Dinput,Iinput)
P1=X(1,1);
P2=X(2,1);
Isc=X(3,1);
Qlocal=X(4,1);
I=X(5,1);
Q=X(6,1);
G=ppval(Gpp,t);

D=Dinput{1,2}(find(Dinput{1,1}<=t,1,'last')); %glucose input
Ib=Iinput{1,2}(find(Iinput{1,1}<=t,1,'last')); %insulin input

Spp=interp1(val.S{1,1},val.Uen,t,'previous');

dSIdt=zeros(8,1); %differential equations used to solve for SI
dSIdt(1,1)=(-val.d1*P1)+D;
dSIdt(2,1)=(-val.d2*P2)+(val.d1*P1);
dSIdt(3,1)=(-val.k2*Isc)+Ib;
dSIdt(4,1)=(-val.k3*Qlocal)+(val.k2*Isc)-(val.kdi*Qlocal);
dSIdt(5,1)=(-I*val.nK)-(val.nL*((I)/(1+(I*val.alpha_I))))-((val.nI/val.Vp)*(I-Q))+(((val.k3*Qlocal)+((1-val.xL)*Spp))/(val.Vp));
dSIdt(6,1)=((val.nI/val.Vq)*(I-Q))-(Q*val.nC);

% dSIdt(7,1)=(G*Q)/(1+(val.alpha_G*Q));
% dSIdt(8,1)=(-val.pg*(G-Greal(1)))+((P2*val.d2)+val.EGP-val.CNS)/(val.Vg));

dSIdt(7,1)=(G*Q)-(Greal(1)*(Ireal(1)/2));
dSIdt(8,1)=(-val.pg*(G-Greal(1)))+(P2*val.d2/val.Vg);
end
