% Function to fit nL and xL.
function P = FitHepaticClearance(P)

global SC GC

tArray = (0 : P.simDuration)';  % Simulation t array [min]

% Time and data arrays.
tI = minutes(P.I.time - P.I.time(1))';  % t of reading [min]
ppI = griddedInterpolant(tI, P.I.value);

%% Analytical Forward Simulation for Q
Q = zeros(length(tArray),1); %analytical solution for Q

% Consider form of dQ/dt = kQ*Q + kI*I.
kQ = -GC.nC - GC.nI/GC.VQ(P, SC); % Constant term coefficent of Q - easier to use
kI = GC.nI/GC.VQ(P, SC);  % Constant term coefficent of Q - easier to use

Q0 = P.I.value(1)/2; %ADM: same assumption as elsewhere
t0 = tArray(1);
Q(1) = Q0;

for ii=2:length(tArray)
    t = tArray(ii);         % Current time value.
    tSpan = tArray(1:ii);   % Array of time values from t0 to t.
    ISpan = ppI(tSpan);     % Array of I values from t0 to t.
    
    % Analytical solution for the Q equation at time==t.
    % Standard soln. for non-homogenous 1st order ODE (page 17).
    Q(ii) = Q0*exp(kQ*(t-t0)) ...
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
% 
%  
% num=length(tArray);
% Qlocal=zeros(num,1); %obtain Qlocal values at each minute from ODE solver solution
% for ii=1:num-1
%     for j=1:length(tFinal)
%         if tFinal(j)>=tArray(ii)
%             Qlocal(ii)=YFinal(j,2);
%             break
%         end
%     end
% end
% Qlocal(num)=YFinal(end,2);

%% Parameter ID of I Equation to find nL/xL (pg. 16)

I = ppI(tArray); % I values over all time.
I0 = ppI(t0);
Uen = P.Uen.value;       % Uen values over all time.

% Set coefficients for MLR.
% Consider dI/dt = kI*I + c1*nL + kIQ*(I-Q) + c2*xL + k
kI = -GC.nK;
kIQ = -GC.nI./GC.VQ(P, SC);
k = Uen/GC.VI(P);

% Therefore, integrating:
% I(t) - I(t0) = kI*int{I} + int{c1}*nL + kIQ*int{I-Q} + int{c2}*xL + int{k}
% Renaming cN = int{c1} and cX = int{c2}
% cN*nL + cX*xL = I(t) - I(t0) - kI*int{I} - kIQ*int{I-Q} - int{k}
cN = cumtrapz(tArray, ...
              -I./(1 + GC.alphaI*I));
cX = cumtrapz(tArray, ...
              -Uen/GC.VI(P));
          
% Assembling MLR system:
% [cN(t), cX(t)] * (nL; xL) = [b(t)]
A = [cN cX];
b = I - I0 ...
      - kI * cumtrapz(tArray, GC.nK*I) ...
      - kIQ * cumtrapz(tArray, I-Q) ...
      - cumtrapz(tArray, k); 

% Solve and save.
x = A\b;
nL = x(1);
xL = x(2);

lb = 1e-7; %"lower bound"
nL = max(nL, lb);
xL = max(xL, lb);

P.nL = nL;
P.xL = xL;
end