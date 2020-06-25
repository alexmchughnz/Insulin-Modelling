function F = MakeDebugPlot(P, DP)
persistent patient;
    if isempty(patient)
        patient = struct();
    end
persistent plotset;
    if isempty(plotset)
        plotset = struct();
    end
persistent set;
    if isempty(set)
        set = 0;
    end
persistent num;
    if isempty(num)
        num = 1;
    end
    
if ~isequal(patient, P)
   % New patient, reset set and num counters.
   set = 1;
   num = 1;
end
if ~isequal(plotset, DP)
   % New set of plots, increment set counter and reset num.
   set = set + 1;
   num = 1;
end

fignum = 100 * P.patientNum + 10*set + 1*num;
F = figure(fignum);

num = num + 1;
plotset = DP;
patient = P;
end
