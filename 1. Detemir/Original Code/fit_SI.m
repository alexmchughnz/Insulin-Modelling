function [sys] = fit_SI(sys)
% function to fit SI for finding SI

%bring blood test & meal data in from from patient file via sys (already
%rea the data)
default_SI = 10.8e-4;
minutes_per_interval = 360; %360 = 6 hourly, this is a good balance
%dunno why im doing this but hey...
optionsShort = odeset('RelTol',1e-5,'AbsTol',1e-4,'MaxStep',0.2,'InitialStep',0.01);    % for first minute
optionsLong = odeset('RelTol',1e-5,'AbsTol',1e-4,'InitialStep',0.1);    % for subsequent minutes

%find how many intervals we are fitting for
intervals_to_fit = floor(sys.sim_max_t/minutes_per_interval);

%create array of hourly SI values
SI = zeros(1, intervals_to_fit);
temp_SI = zeros(1, intervals_to_fit * minutes_per_interval);

%getting index of first and last simulation times
i_start = find(sys.Data.bg3_time == sys.sim_start_t);
i_end = find(sys.Data.bg3_time == sys.sim_end_t);

relative_bg3_times = minutes(sys.Data.bg3_time(i_start : i_end) - sys.Data.bg3_time(i_start));
relative_ins_times = minutes(sys.Data.PlasmaI_time - sys.Data.PlasmaI_time(1));

%fitting piecewise polynomial for glucose and insulin traces
Gpp = interp1(relative_bg3_times, sys.Data.bg3(i_start : i_end), 'linear', 'pp'); %picks only those SI values within the simulation time
Ipp = interp1(relative_ins_times, sys.Data.PlasmaI, 'linear', 'pp');

sys.results.Gpp = Gpp;
sys.results.Ipp = Ipp;

ODE_initials = [0.001, 0, 0, 9, 0, 0];

%setting up uge array to store all the SI values
sys.GC.SI = ones(1, sys.sim_max_t) * default_SI;

start_time = 1;
stop_time = start_time + minutes_per_interval;
ta = 1;
tb = 361;

for ii = 1 : intervals_to_fit    
    % calculate the hourly SI value
    % t = start_time : start_time + 59;
    
    [TI1, Ints1] = ode45(@FAERIES_integrals, (start_time : start_time+1)', ODE_initials, optionsShort, sys, Gpp, Ipp);
    ODE_initials = Ints1(end, :)';
    [TI2, Ints2] = ode45(@FAERIES_integrals, (start_time+1 : start_time+minutes_per_interval-1)', ODE_initials, optionsLong, sys, Gpp, Ipp);
    TI = [TI1(1, :); TI2];
    Ints = [Ints1(1, :); Ints2];
    
%     [sim_time, integrals] = ode45(@FAERIES_integrals, t', ODE_initials, optionsLong, sys, Gpp, Ipp);
%     TI = sim_time;
%     Ints = integrals;

%     ta = start_time; tb = stop_time;
    
    %whats going on here?
%     A = zeros(1,1); b = zeros(1,1);
        %think this might be wrong....
    %do i have to create a whole 360 length array?? probs...
   
%     A(1,1) = Ints(find(TI<=tb, 1, 'last'),5) - Ints(find(TI<=ta, 1, 'last'),5)
%     b(1,1) = -(ppval(Gpp,tb) - ppval(Gpp,ta)) - (Ints(find(TI<=tb, 1, 'last'),6) - Ints(find(TI<=ta, 1, 'last'),6))
    

    A = Ints(:, 5);
    b = -diff(ppval(Gpp, [ta, tb]))' - Ints(:, 6);

    
%     B = Ints( [diff(Gpp), 0]/dt
    %update the conditions
    
    ODE_initials = Ints(end, :);
    %Solve the linear system
    SI(ii) = A\b;
    
    start_time = start_time + minutes_per_interval;
    temp_SI((ii-1)*minutes_per_interval+1 : ii*minutes_per_interval) = SI(ii);
    ta = ta + 360;
    tb = tb + 360;

end

% start_time = minutes(sys.Data.PlasmaI_time(1) - sys.Data.bg3_time(1));
% sys.GC.SI(start_time : start_time+length(temp_SI) - 1) = temp_SI;

sys.GC.SI(1 : length(temp_SI)) = temp_SI;

% stairs(SI)
% title('6 hour fitting')
% ylabel('SI')
% xlabel('multiples of 6 hours')
% axis([1 intervals_to_fit 0 3e-3])


fprintf('SI fit successfully. SI ( 1e-4 ) = ')
disp(SI*1e4)
end

function [dy_dt] = FAERIES_integrals(t, y, sys, Gpp, Ipp)
%sets up the integrals to be solved by fit_SI model
%initial conditions, y0
%y1 - q_sto1
%y2 - q_sto2
%y3 - q_gut
%y4 - Q
EGP = 0.96; %bypass all that nasty shit ---> will need to change this... eventually...

%how many integrals are we doing - 6
dy_dt = zeros(6,1);

%GI model compartments
q_sto1 = y(1);
q_sto2 = y(2);
q_gut = y(3);

G = ppval(Gpp, t); %plasma glucose - using the linear interpolated values for simplicity
I = ppval(Ipp, t); %plasma insulin - this is used to cut out the whole SC compartment method
Q = y(4);

dy_dt(1) = - sys.GI.k21*q_sto1 + eating(sys, t);% + carbohydrates(sys, t); %rate of change of glucose in stomach solid state
dy_dt(2) = - kempt(t,q_sto1, q_sto2,sys)*q_sto2 + sys.GI.k21*q_sto1; %rate of change of glucose in liquid solid state
dy_dt(3) = - sys.GI.kabs*q_gut + kempt(t,q_sto1,q_sto2,sys)*q_sto2; %rate of change of glucose in gut

dy_dt(4) = sys.GC.nI*(I - Q) - (sys.GC.nc * Q)/(1+sys.GC.alphaG*Q); %rate of change of insulin in the interstitium

Ra = sys.GI.f*sys.GI.kabs*q_gut; %rate glucose absorbed into bloodstream

dy_dt(5) = G*Q/(1 + sys.GC.alphaG*Q); %splitting up the main glucose ode (this term is normally multiplied by SI - this becomes A matrix)
dy_dt(6) = (Ra + EGP + infusion(sys, t) - sys.GC.CNS)/sys.GC.VG - sys.GC.pg*G; %this is the second term - forms the first part of b matrix.

end