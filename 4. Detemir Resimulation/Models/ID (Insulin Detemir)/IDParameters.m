% Initial Values
ID.ISC0 = 0;
ID.QDFLocal0 = 0;
ID.QDBLocal0 = 0;
ID.IDF0 = 0;
ID.IDB0 = 0;
ID.QDF0 = 0;
ID.QDB0 = 0;

% Parameters
ID.ka   = 0.0078;   % Hexameric Detemir dissociation rate [1/min]
ID.kd1  = 0.96;     % Detemir unbind fraction [1/min]
ID.kd2  = 0.04;     % Detemir bind fraction [1/min]
ID.kdi  = 0.0594;   % Detemir interstitial degradation [1/min]
ID.kb   = 0.0181;   % Unbound I diffusion from local to blood [1/min]
% NOTE: nK doesn't match Klenner paper (nK = 0).
ID.nK   = 0.06;     % Detemir renal clearance [1/min]
ID.nDL  = 0.147;    % Detamir hepatic insulin clearance rate [1/min]
ID.nDI  = 0.06;     % Detemir trans-endothelial diffusion rate [1/min]
ID.nDC  = 0.032;    % Detemir insulin degredation rate [1/min]