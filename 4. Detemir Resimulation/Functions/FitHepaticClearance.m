% Function to fit nL and xL.
function P = FitHepaticClearance(P)

global SC GC

tArray = (0 : P.simDuration)';  % Simulation t array [min]

% Time and data arrays.
tI = minutes(P.I.time - P.I.time(1))';  % t of reading [min]
ppI = griddedInterpolant(tI, P.I.value);

%% Analytical Forward Simulation for Q
QSolution = zeros(length(tArray),1); %analytical solution for Q

% Consider form of dQ/dt = kQ*Q + kI*I.
kQ = -GC.nC - GC.nI/GC.VQ(P, SC); % Constant term coefficent of Q - easier to use
kI = GC.nI/GC.VQ(P, SC);  % Constant term coefficent of Q - easier to use

Q0 = P.I.value(1)/2; %ADM: same assumption as elsewhere
t0 = tArray(1);
QSolution(1) = Q0;

for ii=2:length(tArray)
    t = tArray(ii);         % Current time value.
    tSpan = tArray(1:ii);   % Array of time values from t0 to t.
    ISpan = ppI(tSpan);     % Array of I values from t0 to t.
    
    % Analytical solution for the Q equation at time==t.
    % Standard soln. for non-homogenous 1st order ODE (page 17).
    QSolution(ii) = Q0*exp(kQ*(t-t0)) ...
                    + trapz(tSpan, kI*exp(kQ*(t - tSpan).*ISpan));
end

%% Incremental solve for Insulin injection parameters over 1 min intervals
% % ADM: I don't think we need to do this for ID, as the I equation does not
% % depend on any of the ID compartments!
% Y0 = [0 0]; %initial conditions of ODEs
% % ADM: this will need to be ID model parameters, won't it?
% 
% options=odeset('RelTol',1e-8,'AbsTol',1e-8);
% 
% 
% tFinal=[];
% YFinal=[];
% for ii = 1 : length(tArray)-1 %loop through ts
%     % Range to solve over for this iteration.
%     tInterval = [tArray(ii) tArray(ii+1)];
%     [t,Y] = ode45(@nLxLModelODE, tInterval, Y0, options, P);
%     Y0 = Y(end, :)'; %sets up initial conditions for next t
%     
%     % Append next minute of results.
%     tFinal = [tFinal(1:end-1);
%               t];
%     YFinal = [YFinal(1:end-1, :);
%               Y];
% end


 %%
 
num=length(tArray);
Qlocal=zeros(num,1); %obtain Qlocal values at each minute from ODE solver solution
for ii=1:num-1
    for j=1:length(tFinal)
        if tFinal(j)>=tArray(ii)
            Qlocal(ii)=YFinal(j,2);
            break
        end
    end
end
Qlocal(num)=YFinal(end,2);

IInterp = ppI(tArray);
Uen = P.Uen.value;
cn=cumtrapz(tArray, IInterp./(1+(GC.alphaI .* IInterp)));
cx=cumtrapz(tArray, Uen/GC.VI(P));
b = (ppI(0) - IInterp) - cumtrapz(tArray, GC.nK*IInterp) ...
      - cumtrapz(tArray,(GC.nI./GC.VI(P)).*(IInterp-QSolution)) ...
      + cumtrapz(tArray,(SC.k3./GC.VI(P)).*Qlocal) ...
      + cumtrapz(tArray, Uen/GC.VI(P));
% NOTE: SC.k3 almost definitely wrong here, need to figure out correct value
% in conjunction with new SC model.
  
  A=[cn cx];

x=A\b;

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