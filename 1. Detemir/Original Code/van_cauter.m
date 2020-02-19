function[Uen] = van_cauter(sys)


%van cauter model
 
%input system struct
%output Uen
 
%here we go...
 
t = minutes(sys.Data.Cpep_time - sys.Data.Cpep_time(1));
C = sys.Data.Cpep;
 
%hard coded pt 1 data
% t = [0, 115, 160, 235, 290, 315, 370, 580, 640, 670, 730, 970, 1425, 1495, 1525, 1590, 1690, 1725, 1760, 1810, 2005, 2080, 2110, 2180];
% C = [3700, 4360, 3310, 3060, 5340, 4930, 3930, 3610, 4400, 5220, 5580, 1800, 1440, 3040, 4970, 5780, 3430, 4370, 5910, 5720, 2280, 3470, 4840, 4340]; %concentration in pmol/L
 
%initialising rate constants
k1 = 0.0478;    %+/- 0.019
k2 = 0.0516;    %+/- 0.013

%volume of distribution
V = 5; %+/- 0.83
 
C = C * V; %converting units from pmol/L to pmol
C = interp1(t, C, 'linear', 'pp');
 
%making vector of 1min spacing
dt = 1;
t = 0 : dt : t(end);
n = length(t);
C = ppval(C, t);

%making time dependant clearance rate
k3 = zeros(1, n);
k3(1 : sys.GC.t_nk_change) = sys.GC.nK/2;
k3(sys.GC.t_nk_change + 1 : end) = sys.GC.nK/2;
k1 = k1 * ones(1, n);

%initialising Y and S vectors
Y = zeros(1, n);
 
%getting first Y value - assumes dY/dt = 0
Y(1) = k1(1)/k2*C(1);
Y = exp(-k2*t).*(Y(1) + cumtrapz(t,exp(k2*t).*k1.*C));


%calculating endogenous secretion of insulin --- this is a rate
S = ( [diff(C), 0]/dt + (k1 + k3).*C - k2*Y);
% Uen2 = ([diff(C), 0]/dt + (k1 + k3)*C - k2*Y2) * V; % p. docherty's one

Uen = S * 1e-12 * 5808 / 33.7e-6 * 1000; %convert from pmol/min to mU/min

% C = C/V;
% Y = Y/V;
% Uen = Uen * 1e-12 * 5808 / 28.8; %convert to U/L/min
% total_secretion = cumtrapz(Uen/1000);

% figure(1)
% plot(t, C)
% hold on
% plot(t, Y)
% title('C-peptide in the periphery and central compartments')
% xlabel('time, min')
% ylabel('C-peptide, pmol/L')
% legend('Central', 'Periphery')
% hold off
%
% figure(2)
% hold on
% plot(t, Uen)
% title('Endogenous secretion of insulin')
% xlabel('time, min')
% ylabel('Endogenous insulin secretion rate, mU/min')
% umax = 267;
% plot([0 t(end)], [umax umax])
% 
% figure(3)
% plot(t, S)
% title('Endogenous secretion of insulin')
% xlabel('time, min')
% ylabel('Endogenous insulin secretion rate, pmol/min')
% 
% figure(4)
% plot(total_secretion)
% title('Total production of insulin')
% xlabel('time, min')
% ylabel('Endogenous insulin, U')

disp('Fit Uen from C-peptide data')
end