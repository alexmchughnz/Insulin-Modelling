function F = MakeDebugPlot(figTitle, P, DP)
% Generates a debug plot object.
% INPUTS:
%   figTitle - title of figure
%   P        - patient struct
%   DP       - debug plot group
% OUTPUT:
%   F  - figure handle for debug plot

global DEBUGPLOTS

GetFigNum = @(p, s, n) 100*p + 10*s + n;


persistent seenPatients;
if isempty(seenPatients)
    seenPatients = [];
end

persistent stateData;
if isempty(stateData)
    stateData = {};
end

    
if ~ismember(P.patientNum, seenPatients)
   % New patient.
   seenPatients = [seenPatients P.patientNum];
   set = 1;
   num = 1;
elseif ~isequal(stateData{P.patientNum}.prevDP, DP)    
   % New set of plots for patient in previous set.
   set = stateData{P.patientNum}.prevset + 1;
   num = 1;
else
   % New plot in previous set.
   set = stateData{P.patientNum}.prevset;
   num = stateData{P.patientNum}.prevNum + 1;
end

% Make figure.
fignum = GetFigNum(P.patientNum, set, num);
F = figure(fignum);
axes
hold on
title(F.Children(1), P.patientCode + ": " + figTitle);
DEBUGPLOTS.FIGURES = vertcat(DEBUGPLOTS.FIGURES, [P.patientNum, set, num]);

% Update pointers.
stateData{P.patientNum}.prevDP = DP;
stateData{P.patientNum}.prevset = set;
stateData{P.patientNum}.prevNum = num;
end
