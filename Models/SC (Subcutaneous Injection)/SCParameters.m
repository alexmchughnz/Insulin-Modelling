function SC = SCParameters(~)

% Initial Values
SC.ISC0 = 0;      % Subcutaneous insulin [mU]
SC.QLocal0 = 0;   % Local interstitium insulin [mU]

% Parameters
SC.kdi = 0.006;   % Breakdown of insulin in interstitium [1/min]
SC.ks2  = 0.0104;  % Diffusion from SC space to interstitium [1/min]
SC.ks3  = 0.060;   % Absorption of insulin into plasma [1/min]

end