function [] = PlotResults(results)
% Plots patient data:
%     - blood glucose, G
%     - plasma insulin, I
%     - insulin sensitivity, SI
%     - estimated endogenous insulin secretion, Uen
% INPUTS:
%   P - patient struct


%% Glucose
figure(1)
plot(results.tArray, results.G);

title('Plasma Glucose')
xlabel('Time [min]')
ylabel('Plasma Glucose, G [mmol/L]')

end

