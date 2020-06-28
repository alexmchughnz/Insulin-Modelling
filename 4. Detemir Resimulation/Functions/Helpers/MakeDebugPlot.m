function F = MakeDebugPlot(P, DP)

global DEBUGPLOTS

GetFigNum = @(p, s, n) 100*p + 10*s + n;

persistent patient;
    if isempty(patient)
        patient = struct('patientNum', 0);
    end
persistent plotset;
    if isempty(plotset)
        plotset = struct();
    end
persistent set;
    if isempty(set)
        set = 1;
    end
persistent num;
    if isempty(num)
        num = 1;
    end
    
if ~isequal(patient.patientNum, P.patientNum)
   % New patient, reset set and num counters.
   set = 1;
   num = 1;
elseif ~isequal(plotset, DP)
   % New set of plots, increment set counter and reset num.
   set = set + 1;
   num = 1;
end

% Make figure.
fignum = GetFigNum(P.patientNum, set, num);
F = figure(fignum);
DEBUGPLOTS.FIGURES = vertcat(DEBUGPLOTS.FIGURES, [P.patientNum, set, num]);

% Update pointers.
num = num + 1;
plotset = DP;
patient = P;
end
