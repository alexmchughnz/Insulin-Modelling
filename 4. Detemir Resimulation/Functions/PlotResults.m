function [] = PlotResults(P, source)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P      - patient struct
%   source - string of trial type, e.g. "DISST"

global C

ToDateTime = @(mins) P.data.simTime(1) + minutes(mins);

tArray = P.results.tArray;     % Time of results [min]
dtArray = ToDateTime(tArray);  % ''              [datetime]

% Set up figure.
patientLabel = sprintf("Patient %d: ", P.patientNum);
F = PanelFigures(3, 3);


%% Glucose
subplot(4, 1, 1)
hold on

[tG, vG] = GetSimTime(P, P.data.G);
dtG = ToDateTime(tG);
plt = plot(dtG, vG, 'r*');
plt.DisplayName = 'Blood Test';

ppG = griddedInterpolant(tG, vG);
plt = plot(dtArray, ppG(tArray), 'b');
plt.LineWidth = 1;
plt.DisplayName = 'Interpolation';

plt = plot(dtArray, P.results.G, 'k');
plt.DisplayName = 'Model Prediction';

lineBounds = ylim;
for ii = 1:length(P.results.nLxLFitBounds)
    split = ToDateTime(P.results.nLxLFitBounds(ii));
    L = line([split split], lineBounds);
    L.LineWidth = 0.5;
    L.Color = 'k';
    L.HandleVisibility = 'off';
end

title([patientLabel 'Plasma Glucose'])
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')
legend()

if source == "Detemir"
    datetick('x')
end
ylim([4 15])


%% Glucose Error
ax = subplot(4, 1, 2);

iiG = GetTimeIndex(tG, tArray);
simG = P.results.G(iiG);
GError = 100*abs((simG - vG) ./ vG);
plot(dtG, GError, 'r');

title([patientLabel 'Plasma Glucose Error'])
xlabel('Time')
ylabel('Error [\%]')


%% Insulin (+ Detemir)
subplot(4, 1, 3)
hold on
if source == "Detemir"
    [tI, vI] = GetSimTime(P, P.data.ITotal);  % [pmol/L]
    tI =ToDateTime(tI);
    
    ppI = griddedInterpolant(tI, vI);
    
    I = C.mU2pmol(P.results.I + P.results.IDF);  % [mU/L] -> [pmol/L]
    
    plttitle = [patientLabel 'Plasma Insulin + Detemir'];
    pltxlabel = 'Time';
    pltylabel = 'Plasma Insulins, I + IDF [pmol/L]';
    pltxarray = dtArray;
    
    datetick('x')
    
elseif source == "DISST"
    [tI, vI] = GetSimTime(P, P.data.I);  % [pmol/L]
    
    ppI = griddedInterpolant(tI, vI);
    
    I = C.mU2pmol(P.results.I);  % [mU/L] -> [pmol/L]
    
    plttitle = [patientLabel 'Plasma Insulin'];
    pltxlabel = 'Time';
    pltylabel = 'Plasma Insulin, I [pmol/L]';
    pltxarray = tArray;
end

lineBounds = ylim;
for ii = 1:length(P.results.nLxLFitBounds)
    split = ToDateTime(P.results.nLxLFitBounds(ii));
    L = line([split split], lineBounds);
    L.LineWidth = 0.5;
    L.Color = 'k';
    L.HandleVisibility = 'off';
end

plt = plot(tI, vI, 'r*');
plt.DisplayName = 'Blood Test';

plt = plot(pltxarray, ppI(tArray), 'b');
plt.LineWidth = 1;
plt.DisplayName = 'Interpolation';

plt = plot(pltxarray, I, 'k');
plt.DisplayName = 'Model Prediction';


title(plttitle)
xlabel(pltxlabel)
ylabel(pltylabel)
legend()

%% Insulin Error
ax = subplot(4, 1, 4);

iiI = GetTimeIndex(tI, tArray);
simI = I(iiI);
ITotalError = 100*abs((simI - vI) ./ vI);
plot(tI, ITotalError, 'r');

title([patientLabel 'Plasma Insulins Error'])
xlabel('Time')
ylabel('Error [\%]')


%%
path = fullfile("Plots", "patient" + P.patientNum);
savefig(F, path);
fprintf("P%d: Plotted results.\n", P.patientNum)

end

