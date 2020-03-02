%% function that solves for  G, I, and Q
function [output,T_GIQ]=GIQ_Model(val,Greal,Ireal,timestep,Dinput,Iinput)
format long

ic=zeros(1,7); %initial conditions of ODEs
ic(1,1)=0; %P1
ic(1,2)=0; %P2
ic(1,3)=0; %Isc
ic(1,4)=0; %Qlocal
ic(1,5)=Greal(1); %G
ic(1,6)=Ireal(1); %I
ic(1,7)=Ireal(1)/2; %Q

options=odeset('RelTol',1e-8,'AbsTol',1e-8);
tfinal=[];
icfinal=[];
for i=1:length(timestep)-1 %loop through timesteps
    [t,ysol]=ode45(@GIQ_Model_ODE,[timestep(i) timestep(i+1)],ic,options,val,Dinput,Iinput,Greal,Ireal);
    ic=ysol(end,:)'; %sets up initial conditions for next timestep
    tfinal=[tfinal(1:end-1,:);t]; %put results onto the end of the last timestep
    icfinal=[icfinal(1:end-1,:);ysol];
end
output=icfinal;
T_GIQ=tfinal;
end

%% function that sets up differential equations to solve for G, I, and Q
function dGIQdt=GIQ_Model_ODE(t,X,val,Dinput,Iinput,Greal,Ireal)
P1=X(1,1);
P2=X(2,1);
Isc=X(3,1);
Qlocal=X(4,1);
G=X(5,1);
I=X(6,1);
Q=X(7,1);

D=Dinput{1,2}(find(Dinput{1,1}<=t,1,'last')); %glucose input
Ib=Iinput{1,2}(find(Iinput{1,1}<=t,1,'last')); %insulin input

Spp=interp1(val.S{1,1},val.Uen,t,'previous');

dGIQdt=zeros(7,1); %differential equations used to solve for GIQ
dGIQdt(1,1)=(-val.d1*P1)+D; %P1
dGIQdt(2,1)=(-val.d2*P2)+(val.d1*P1); %P2
dGIQdt(3,1)=(-val.k2*Isc)+Ib; %Isc
dGIQdt(4,1)=(-val.k3*Qlocal)+(val.k2*Isc)-(val.kdi*Qlocal); %Qlocal
% dGIQdt(5,1)=(-val.pg*(G-Greal(1)))-((val.SI*G*Q)/(1+(Q*val.alpha_G)))+(((val.d2*P2)+val.EGP-val.CNS)/(val.Vg)); %G
dGIQdt(5,1)=(-val.pg*(G-Greal(1)))-(val.SI*((G*Q)-(Greal(1)*(Ireal(1)/2))))+(P2*val.d2/val.Vg); %G
dGIQdt(6,1)=(-I*val.nK)-(val.nL*((I)/(1+(I*val.alpha_I))))-((val.nI/val.Vp)*(I-Q))+(((val.k3*Qlocal)+((1-val.xL)*Spp))/(val.Vp)); %I
dGIQdt(7,1)=((val.nI/val.Vq)*(I-Q))-(Q*val.nC); %Q

end