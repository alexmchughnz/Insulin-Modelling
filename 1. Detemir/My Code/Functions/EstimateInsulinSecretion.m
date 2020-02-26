function P = EstimateInsulinSecretion(GC, P)
% Estimates pancreatic insulin secretion rate (Uen) using a patient's 
% C-peptide data.
% Model from Van Cauter et al. (1992).
% Code based on van_cauter.m.
% INPUTS:
%   GC  - glycaemic control parameter set
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen

% STUB: Use previous Van Cauter function

    %van cauter model

    %input system struct
    %output Uen

    %here we go...

    t = minutes(P.CPep.time - P.CPep.time(1));
    C = P.CPep.value;

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
    t = (0 : dt : t(end))';
    n = length(t);
    C = ppval(C, t);

    %making time dependant clearance rate
    k3 = zeros(n, 1);
    k3(1 : GC.nKChangeTime) = GC.nK/2;
    k3(GC.nKChangeTime + 1 : end) = GC.nK/2;
    k1 = k1 * ones(n,1);

    %initialising Y and S vectors
    Y = zeros(n,1);

    %getting first Y value - assumes dY/dt = 0
    Y(1) = k1(1)/k2*C(1);
    Y = exp(-k2*t).*(Y(1) + cumtrapz(t,exp(k2*t).*k1.*C));


    %calculating endogenous secretion of insulin --- this is a rate
    S = ( [diff(C); 0]/dt + (k1 + k3).*C - k2*Y);
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

    fprintf('P%d: Estimated Uen.\n', P.patientNum);    
%\STUB

% Write value to patient struct.
P.Uen.value = Uen;
P.Uen.time  = t;



% % Extract C-peptide data.
% C = P.CPep.value;  % [pmol/L]
% t = P.CPep.time;   % [datetime]
% 
% % Tidy data and change units.
% C = C * GC.VI(P);       % Amount of C-peptide [pmol]
% t = minutes(t - t(1));  % Time from first reading [min]
% 
% 
% % Interpolate data over 1 min spacings.
% C = griddedInterpolant(t, C);
% 
% 
% t = 0 : t(end);         %   -spaced in 1 min grid
% 
% 
% 
% 
% 
% 
% 
% % Rate constants.
% % ***Currently taken from van_cauter.m***
% 
% k1 = 0.0478;    %+/- 0.019
% k2 = 0.0516;    %+/- 0.013
% k3 = GC.nK;     % NEEDS TO BE HALVED AT SOME POINT?



end

