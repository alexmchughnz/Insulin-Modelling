function [] = PlotResults(variants)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P - patient struct


%% Glucose
figure(1)
hold on
for ii = 1:length(variants)
    V = variants{ii};
    results = V.results;
    
    p = plot(results.tArray, results.G);
    p.DisplayName = "SI = " + V.SI;
end

legend

title('Plasma Glucose')
xlabel('Time [min]')
ylabel('Plasma Glucose, G [mmol/L]')

end

