%% function that fits nL and xL
function [nL,xL]=nLxL_Model(Treal,Ireal,val,timestep,Iinput)
format long

time=[timestep(1):1:timestep(end)]'; %create a time vector with 1 min intervals

Ipp=interp1(Treal,Ireal,time); %fit a piecewise fxn to the insulin measurements
K=val.nC+(val.nI/val.Vq); %constant term - easier to use
Qpp=zeros(length(time),1); %analytical solution for Q
Q0=Ireal(1)/2;
Qpp(1,1)=Q0;
for i=2:length(time)
    Qpp(i,1)=Q0.*exp(-K.*(time(i)+abs(time(1))))+trapz(time(1:i,1),(val.nI./val.Vq).*Ipp(1:i,1).*exp(-K.*(time(i,1)-time(1:i,1))));
end

ic=zeros(1,2); %initial conditions of ODEs
ic(1,1)=0; %Isc
ic(1,2)=0; %Qlocal

options=odeset('RelTol',1e-8,'AbsTol',1e-8);

tfinal=[];
icfinal=[];
for i=1:length(timestep)-1 %loop through timesteps
    [t,ysol]=ode45(@nLxL_Model_ODE,[timestep(i) timestep(i+1)],ic,options,val,Iinput);
    ic=ysol(end,:)'; %sets up initial conditions for next timestep
    tfinal=[tfinal(1:end-1,:);t]; %put results onto the end of the last timestep
    icfinal=[icfinal(1:end-1,:);ysol];
end

num=length(time);
Qlocal=zeros(num,1); %obtain Qlocal values at each minute from ODE solver solution
for i=1:num-1
    for j=1:length(tfinal)
        if tfinal(j)>=time(i)
            Qlocal(i)=icfinal(j,2);
            break
        end
    end
end
Qlocal(num)=icfinal(end,2);

cn=cumtrapz(time,Ipp./(1+(val.alpha_I.*Ipp)));
cx=cumtrapz(time,val.Uen./val.Vp);
c=(Ipp(1,1)-Ipp)-cumtrapz(time,val.nK.*Ipp)-cumtrapz(time,(val.nI./val.Vp).*(Ipp-Qpp))+cumtrapz(time,(val.k3./val.Vp).*Qlocal)+cumtrapz(time,val.Uen./val.Vp);
A=[cn cx];
x=A\c;
nL=x(1,1);
xL=x(2,1);

lb=1e-7;
nL=max(nL,lb);
xL=max(xL,lb);
end

%% function that sets up differential equations to solve for Isc and Qlocal
function dnLxLdt=nLxL_Model_ODE(t,X,val,Iinput)
Isc=X(1,1);
Qlocal=X(2,1);

Ib=Iinput{1,2}(find(Iinput{1,1}<=t,1,'last')); %insulin input

dnLxLdt=zeros(2,1); %differential equations used to solve for GIQ
dnLxLdt(1)=(-val.k2*Isc)+Ib; %Isc
dnLxLdt(2)=(-val.k3*Qlocal)+(val.k2*Isc)-(val.kdi*Qlocal); %Qlocal
end