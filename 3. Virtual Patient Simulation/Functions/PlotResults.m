function [] = PlotResults(variants)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen


%% Glucose
figure(1)
hold on

% Hypo/hyperglycaemic limits.
rectangle('Position', [0 0 120 3.5], ...
          'FaceColor', '#d58a94', ...
          'LineWidth', 1e-3)
rectangle('Position', [0 3.5 120 0.5], ...
          'FaceColor', '#ffd1df', ...
          'LineWidth', 1e-3)
line([0 150], [8 8], ...
     'Color', 'k', ...
     'LineWidth', 0.01, ...
     'HandleVisibility', 'off')
      
% Plot data.
for ii = 1:length(variants)
    V = variants{ii};
    results = V.results;
    
    p = plot(results.tArray, results.G);
    p.DisplayName = "SI = " + V.SI;
end


ylim([0 15])
legend
grid on

title('Plasma Glucose')
xlabel('Time [min]')
ylabel('Plasma Glucose, G [mmol/L]')

end

