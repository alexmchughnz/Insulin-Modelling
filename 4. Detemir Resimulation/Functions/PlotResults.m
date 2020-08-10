function [] = PlotResults(P)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P      - patient struct

patientfigure = @(n) figure(10*(10*P.patientNum + n));

tArray = P.results.tArray;     % Time of results [min]

% Set up figure.
patientLabel = sprintf("%s: ", P.patientCode);


%% Plasma Glucose
patientfigure(1)
hold on

[tG, vG] = GetSimTime(P, P.data.G);
plt = plot(tG, vG, 'r.');
plt.DisplayName = 'Blood Test';

ppG = griddedInterpolant(tG, vG);
plt = plot(tArray, ppG(tArray), 'b');
plt.LineWidth = 1;
plt.DisplayName = 'Interpolation';

plt = plot(tArray, P.results.G, 'k');
plt.DisplayName = 'Model Prediction';

title(patientLabel + "Plasma Glucose")
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')
legend()

ylim([4 15])


%% Insulin (+ Detemir)
patientfigure(2)
hold on

if P.source == "Detemir"
    [tI, vI] = GetSimTime(P, P.data.ITotal);  % [mU/L]    
    I = P.results.I + P.results.IDF;  % [mU/L]
    
    plttitle = patientLabel + "Plasma Insulin + Detemir";
    pltylabel = 'Plasma Insulins, I + IDF [mU/L]';
    
    datetick('x')
    
else
    [tI, vI] = GetSimTime(P, P.data.I);  % [mU/L]
    I = P.results.I;  % [mU/L]
    
    plttitle = patientLabel + "Plasma Insulin";
    pltylabel = 'Plasma Insulin, I [mU/L]';
end

plt = plot(tI, vI, 'r.');
plt.DisplayName = 'Blood Test';

% ppI = griddedInterpolant(tI, vI);
% plt = plot(pltxarray, ppI(tArray), 'b');
% plt.LineWidth = 1;
% plt.DisplayName = 'Interpolation';

plt = plot(tArray, I, 'k');
plt.DisplayName = 'Model Prediction';


title(plttitle)
xlabel('Time [min]')
ylabel(pltylabel)
legend()



%%
fprintf("P%d: Plotted results.\n", P.patientNum)

end

