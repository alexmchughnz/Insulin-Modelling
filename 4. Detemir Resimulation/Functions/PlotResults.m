function [] = PlotResults(P)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P      - patient struct

global C

patientfigure = @(n) figure(10*P.patientNum + n);

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

if isfield(P.results, 'nLxLFitBounds')
    lineBounds = ylim;
    for ii = 1:length(P.results.nLxLFitBounds)
        split = tArray(P.results.nLxLFitBounds(ii));
        L = line([split split], lineBounds);
        L.LineWidth = 0.5;
        L.Color = 'k';
        L.HandleVisibility = 'off';
    end
end

title(patientLabel + "Plasma Glucose")
xlabel('Time')
ylabel('Plasma Glucose, G [mmol/L]')
legend()

ylim([4 15])


%% Insulin (+ Detemir)
patientfigure(2)
hold on

if P.source == "Detemir"
    [tI, vI] = GetSimTime(P, P.data.ITotal);  % [pmol/L]
    
    ppI = griddedInterpolant(tI, vI);
    
    I = C.mU2pmol(P.results.I + P.results.IDF);  % [mU/L] -> [pmol/L]
    
    plttitle = patientLabel + "Plasma Insulin + Detemir";
    pltylabel = 'Plasma Insulins, I + IDF [pmol/L]';
    
    datetick('x')
    
elseif P.source == "DISST"
    [tI, vI] = GetSimTime(P, P.data.I);  % [pmol/L]
    
    ppI = griddedInterpolant(tI, vI);
    
    I = C.mU2pmol(P.results.I);  % [mU/L] -> [pmol/L]
    
    plttitle = patientLabel + "Plasma Insulin";
    pltylabel = 'Plasma Insulin, I [pmol/L]';
end

if isfield(P.results, 'nLxLFitBounds')
    lineBounds = ylim;
    for ii = 1:length(P.results.nLxLFitBounds)
        split = tArray(P.results.nLxLFitBounds(ii));
        L = line([split split], lineBounds);
        L.LineWidth = 0.5;
        L.Color = 'k';
        L.HandleVisibility = 'off';
    end
end

plt = plot(tI, vI, 'r.');
plt.DisplayName = 'Blood Test';

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

