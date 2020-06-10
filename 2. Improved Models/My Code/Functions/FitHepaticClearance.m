%% function that fits nL and xL
function P = FitHepaticClearance(P)

global SC GC

tSpan = (0 : P.simDuration-1)';  % Simulation t array [min]

% t and data arrays.
tI = minutes(P.I.time - P.I.time(1))';  % t of reading [min]
ppI = griddedInterpolant(tI, P.I.value);

K = GC.nC + GC.nI/GC.VQ(P); %constant term - easier to use
QAna = zeros(length(tSpan),1); %analytical solution for Q
Q0 = P.I.value(1)/2; %ADM: same assumption as elsewhere
QAna(1) = Q0;
for ii=2:length(tSpan)
    QAna(ii)= Q0.*exp(-K.*(tSpan(ii)+abs(tSpan(1)))) ...
             + trapz(tSpan(1:ii), GC.nI./GC.VQ(P)*ppI(1:ii)'.*exp(-K*(tSpan(ii)-tSpan(1:ii))));
end

ic=zeros(1,2); %initial conditions of ODEs
ic(1,1)=0; %Isc
ic(1,2)=0; %Qlocal

options=odeset('RelTol',1e-8,'AbsTol',1e-8);

tfinal=[];
icfinal=[];
for i=1:length(tSpan)-1 %loop through ts
    tInt = [tSpan(i) tSpan(i+1)];
    [t,ysol]=ode45(@nLxLModelODE, tInt, ic, options, P);
    ic=ysol(end,:)'; %sets up initial conditions for next t
    tfinal=[tfinal(1:end-1,:); t]; %put results onto the end of the last t
    icfinal=[icfinal(1:end-1,:); ysol];
end

num=length(tSpan);
Qlocal=zeros(num,1); %obtain Qlocal values at each minute from ODE solver solution
for i=1:num-1
    for j=1:length(tfinal)
        if tfinal(j)>=tSpan(i)
            Qlocal(i)=icfinal(j,2);
            break
        end
    end
end
Qlocal(num)=icfinal(end,2);

IInterp = ppI(tSpan);
Uen = P.Uen.value;
cn=cumtrapz(tSpan, IInterp./(1+(GC.alphaI .* IInterp)));
cx=cumtrapz(tSpan, Uen/GC.VI(P));
c = (ppI(0) - IInterp) - cumtrapz(tSpan, GC.nK*IInterp) ...
      - cumtrapz(tSpan,(GC.nI./GC.VI(P)).*(IInterp-QAna)) ...
      + cumtrapz(tSpan,(SC.k3./GC.VI(P)).*Qlocal) + cumtrapz(tSpan, Uen/GC.VI(P));
% NOTE: SC.k3 almost definitely wrong here, need to figure out correct value
% in conjunction with new SC model.
  
  A=[cn cx];

x=A\c;

nL = x(1);
xL = x(2);

lb=1e-7; %"lower bound"
nL=max(nL,lb);
xL=max(xL,lb);

P.nL = nL;
P.xL = xL;
end

%% function that sets up differential equations to solve for Isc and Qlocal
function dY = nLxLModelODE(t, Y, P)

global SC 

Isc=Y(1,1);
Qlocal=Y(2,1);

IB = P.IBolus(t);

dY=zeros(2,1); %differential equations used to solve for GIQ
dY(1)=(-SC.k2*Isc)+IB; %Isc
dY(2)=(-SC.k3*Qlocal)+(SC.k2*Isc)-(SC.k3*Qlocal); %Qlocal
end