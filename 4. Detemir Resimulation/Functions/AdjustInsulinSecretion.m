function P = AdjustInsulinSecretion(P, method)
% Adjusts pancreatic insulin secretion rate (Uen) within reasonable range
% to provide a more dynamic profile.
% Ranges based on Ormsbee et al. (2020).
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen


global DEBUGPLOTS

%% Setup
variation = 15/100;

tArray = P.results.tArray;  
Uen = P.results.Uen;


%% Generate New Uen Profile    
if method == "alternate"
    factors = repmat([1-variation; 1+variation], ceil(length(Uen)/2), 1);
    factors = factors(1:length(Uen));
    
    newUen = P.results.Uen .* factors;    
end
    

AUC = cumtrapz(tArray, Uen);
newAUC = cumtrapz(tArray, newUen);
    
    
%% Debug Plots
DP = DEBUGPLOTS.AdjustInsulinSecretion;
if DP.Uen
   MakeDebugPlot(P, DP);
   hold on
   
   plt = area(tArray, Uen, ...
       'EdgeColor', 'b', ...
       'FaceAlpha', 0.2);
   plt.DisplayName = sprintf("Measured $U_{en}$ (AUC = %.3g)", AUC(end));
   
   plt = area(tArray, newUen, ...
       'EdgeColor', 'r', ...
       'FaceAlpha', 0.2);
   plt.DisplayName = sprintf("Adjusted $U_{en}$ (AUC = %.3g)", newAUC(end));
   
   plt = plot(tArray, Uen*(1+variation), 'k--', 'LineWidth', 1);
   plt.HandleVisibility = 'off';   
   plt = plot(tArray, Uen*(1-variation), 'k--', 'LineWidth', 1);
   plt.HandleVisibility = 'off';
   
   title(sprintf("%s: Adjusted Uen", P.patientCode))   
   xlabel("Time [min]")
   ylabel("Uen [mU/min]")
   legend
end

if DP.AUC    
   MakeDebugPlot(P, DP);
   hold on
   
   
   plt = plot(tArray, AUC, 'b');
   plt.DisplayName = sprintf("Measured $U_{en}$ (AUC = %.3g)", AUC(end));
   
   plt = plot(tArray, newAUC, 'r');
   plt.DisplayName = sprintf("Adjusted $U_{en}$ (AUC = %.3g)", newAUC(end));
   
   
   title(sprintf("%s: Uen AUC", P.patientCode))   
   xlabel("Time [min]")
   ylabel("AUC [mU]")
   legend
end

end

