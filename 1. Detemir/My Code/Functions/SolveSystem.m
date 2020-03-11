function P = SolveSystem(P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

global C

G0Indices = (P.G{3}.time == P.simTime(1));  % Start of G measurements [indices]
G0 = P.G{3}.value(G0Indices);
I0 = C.pmol2mIU(P.I.value(1)); % [pmol/L] -> [mIU/L]

t = (0 : P.simDuration-1)';
Y0 = [0.001;  % qSto1(t=0)
      0; 	  % qSto2(t=0)
      0;      % qGut(t=0)
      0;      % ISC(t=0)
      0;      % QDFLocal(t=0)
      0;      % QDBLocal(t=0)
      0;      % IDF(t=0)
      0;      % IDB(t=0)
      0;      % QDF(t=0)
      0;      % QDB(t=0)
      G0;     % G(t=0)
      I0;     % I(t=0)
      0];     % Q(t=0)

options = odeset('RelTol',1e-5, ...
                 'AbsTol',1e-4);
[t, Y] = ode45(@SystemODE, t, Y0, options, P);

% Write out results.
P.results.time = P.simTime(1) + t/24/60;
P.results.qSto1 = Y(:,1);
P.results.qSto2 = Y(:,2);
P.results.qGut = Y(:,3);
P.results.IDH = Y(:,4);
P.results.QDFLocal = Y(:,5);
P.results.QDBLocal = Y(:,6);
P.results.IDF = Y(:,7);  % Why * 18??????????
P.results.IDB = Y(:,8)*18;
P.results.QDF = Y(:,9)*18;
P.results.QDB = Y(:,10)*18;
P.results.G = Y(:,11);
P.results.I = C.mIU2pmol(Y(:,12)); % [pmol]
P.results.Q = Y(:,13);
end


% Adapted from eq_solve/odefun.
function [dY] = SystemODE(t, Y, P)

global GI ID GC SC

% Assign incoming variables.
qSto1     = Y(1);  % GI model.
qSto2     = Y(2);
qGut      = Y(3);

IDH       = Y(4);  % ID model.
QDFLocal  = Y(5);
QDBLocal  = Y(6);
IDF       = Y(7);
IDB       = Y(8);
QDF       = Y(9);
QDB       = Y(10);

G         = Y(11);  % GC model.
I         = Y(12);
Q         = Y(13);

% Retrieve required derived values.
kEmpt = GetStomachEmptyingRate(t, (qSto1+qSto2), P);
D = GetGlucoseDelivery(t, P); % Current glucose delivery rate [g/min]
Ra = GI.f * GI.kAbs * qGut;   % Rate of glucose appearance in plasma

% Get time-dependent values.
GFast = P.GFast{1+(t>1000)};            % HACKY, rewrite!!!
SI = P.SI(1 + floor(t));
Uen = P.Uen.value(1 + floor(t));

% Gastrointestinal model DEs.
dqSto1 = D - GI.k21*qSto1;
dqSto2 = GI.k21*qSto1 - kEmpt*qSto2;
dqGut  = kEmpt*qSto2 - GI.kAbs*qGut;

% Insulin Detemir model DEs.
dIDH = -ID.ka*IDH; %+ IBolus (TODO)
dQDFLocal = ID.ka*IDH - QDFLocal*(ID.kb + ID.kdi) ...
                - ID.C*(ID.kd1*QDFLocal - ID.kd2*QDBLocal);
dQDBLocal = ID.C*(ID.kd1*QDFLocal - ID.kd2*QDBLocal);
dIDF = ID.kb/GC.VI(P)*QDFLocal - IDF*(ID.nDL + ID.nK) ...
           - ID.nDI*(IDF - QDF) - ID.C*(ID.kd1*QDF - ID.kd2*QDB);
dIDB = ID.C*(ID.kd1*IDF - ID.kd2*IDB);
dQDF = -ID.nDC*QDF + ID.nDI*(IDF - QDF) - ID.C*(ID.kd1*QDF - ID.kd2*QDB);
dQDB = ID.C*(ID.kd1*QDF - ID.kd2*QDB);

% Glycaemic control model DEs.
dG = -GC.pg*(G - GFast) - SI*G*(Q+QDF)/(1 + GC.alphaG*(Q+QDF)) ... % NOTE: Removed infusion term.
          + (Ra + GC.EGP - GC.CNS)/GC.VG(P);
dI = -GC.nK(t)*I - GC.nL*I/(1 + GC.alphaI*I) - GC.nI*(I - Q) + Uen*(1 - GC.xL)/GC.VI(P); %rate of change of plasma insulin concentration
dQ = GC.nI*(I - Q) - GC.nC * Q/(1+GC.alphaG*Q);

% Pack up outgoing variables.
dY = [dqSto1;
      dqSto2;
      dqGut;
      dIDH;
      dQDFLocal;
      dQDBLocal;
      dIDF;
      dIDB;    
      dQDF;     
      dQDB;
      dG;
      dI;
      dQ];

end