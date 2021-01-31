function P = AdjustInsulinSecretion(P, method, time)
% Adjusts pancreatic insulin secretion rate (Uen) within reasonable range
% to provide a more dynamic profile.
% Ranges based on Ormsbee et al. (2020).
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with Uen

DP = DebugPlots().AdjustInsulinSecretion;

%% Setup
relError = @(a, b) abs(a-b)/a;

variation = 15/100;

tArray = P.results.tArray;  
Uen = P.results.Uen;


%% Generate New Uen Profile    
if method == "alternate"
    factors = 1 + variation*repmat([-1; +1], ceil(length(Uen)/2), 1);
elseif method == "sawtooth"
    factors = 1 + variation*repmat([-1; 0; +1], ceil(length(Uen)/2), 1);
elseif method == "square"
    factors = 1 + variation*repmat([-1; -1; +1; +1], ceil(length(Uen)/2), 1);
elseif method == "slowalternate"
    factors = 1 + variation*repmat([-1; -0.5; 0; +0.5; +1; +0.5; 0; -0.5], ceil(length(Uen)/2), 1);
else
    iiAdjust = GetTimeIndex(time, P.results.tArray);    
    if method == "onetooth"
        factors = ones(size(Uen));
        
        factors(iiAdjust) = 1 + variation;
        factors(iiAdjust+1) = 1 - variation;
    end
end

factors = factors(1:length(Uen));
newUen = P.results.Uen .* factors;    


%% Adjust points to conserve AUC.
for ii = 2 : length(Uen)    
    AUCCum = trapz(tArray(1:ii), Uen(1:ii));
    newAUCCum = trapz(tArray(1:ii), newUen(1:ii));
    
    while relError(AUCCum, newAUCCum) > 1/100
        % Pick target to move towards.
        if newAUCCum > AUCCum                 % Need to reduce area!
            sgn = -1;
        else                                  % Need to increase area!
            sgn = +1;
        end
        
        % Move towards target.
        jump = sgn * abs(newAUCCum - AUCCum);
        newUen(ii) = newUen(ii) + jump;        
        
        if newUen(ii) > (1+variation)*Uen(ii)
            newUen(ii) = (1+variation)*Uen(ii);
            break
        elseif newUen(ii) < (1-variation)*Uen(ii)
            newUen(ii) = (1-variation)*Uen(ii);
            break
        end
        
        % Recalculate AUCs.
        AUCCum = trapz(tArray(1:ii), Uen(1:ii));
        newAUCCum = trapz(tArray(1:ii), newUen(1:ii));   
    end    
end

AUC = cumtrapz(Uen);
newAUC = cumtrapz(newUen);
AUCTotal = AUC(end);
newAUCTotal = newAUC(end);

%% Save Result
P.results.Uen = newUen;

    
%% Debug Plots
if DP.Uen
   MakeDebugPlot('Adjusted Uen', P, DP);
   
   plt = area(tArray, Uen, ...
       'EdgeColor', 'b', ...
       'FaceAlpha', 0.2);
   plt.DisplayName = sprintf("Measured $U_{en}$ (AUC = %.3g)", AUCTotal);
   
   plt = area(tArray, newUen, ...
       'EdgeColor', 'r', ...
       'FaceAlpha', 0.2);
   plt.DisplayName = sprintf("Adjusted $U_{en}$ (AUC = %.3g)", newAUCTotal);
   
   plt = plot(tArray, (1+variation)*Uen, 'k--', 'LineWidth', 1);
   plt.HandleVisibility = 'off';   
   plt = plot(tArray, (1-variation)*Uen, 'k--', 'LineWidth', 1);
   plt.HandleVisibility = 'off';
   
 
   xlabel("Time [min]")
   ylabel("Uen [mU/min]")
   legend
end

if DP.AUC    
   MakeDebugPlot('Uen AUC', P, DP);
   
   plt = plot(tArray, AUC, 'b');
   plt.DisplayName = sprintf("Measured $U_{en}$ (AUC = %.3g)", AUCTotal);
   
   plt = plot(tArray, newAUC, 'r');
   plt.DisplayName = sprintf("Adjusted $U_{en}$ (AUC = %.3g)", newAUCTotal);
    
   xlabel("Time [min]")
   ylabel("AUC [mU]")
   legend
end

end

