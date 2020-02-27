function P = SolveSystem(GI, ID, GC, SC, P)
% Fits data to find SI over time for a patient.
% INPUTS:
%   GC  - glycaemic control parameter set
%   GI  - gastrointestinal parameter set
%   P   - patient struct
% OUTPUT:
%   P   - modified patient struct with SI

pmol2mu = @(pmol) pmol * 1e-12 * 5808 / 33.7e-6 * 1000;

G0Indices = (P.G{3}.time == P.simTime(1)); % Start of G measurements [indices]
G0 = P.G{3}.value(G0Indices);
I0 = pmol2mu(P.I.value(1));

t = (0 : P.simDuration-1)';
Y0 = [0.001;            % qSto1(t=0)
      0;                % qSto2(t=0)
      0;                % qGut(t=0)
      0;                % ISC(t=0)
      0;                % QDFLocal(t=0)
      0;                % QDBLocal(t=0)
      0;                % IDF(t=0)
      0;                % IDB(t=0)
      0;                % QDF(t=0)
      0;                % QDB(t=0)
      G0;               % G(t=0)
      I0;               % I(t=0)
      0];               % Q(t=0)

optionsLong = odeset('RelTol',1e-5, ...
                     'AbsTol',1e-4);
[t, Y] = ode45(@SystemODE, t, Y0, optionsLong, GI, ID, GC, SC, P);

P.results = [t Y];

end


% Adapted from eq_solve/odefun.
function [dY] = SystemODE(t, Y, GI, ID, GC, SC, P)

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
kEmpt = GetStomachEmptyingRate(t, (qSto1+qSto2), GI, P);
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
dISC = SC.k2*IDH; %+ UTOTAL - TODO
dQDFLocal = SC.k2*IDH - QDFLocal*(SC.k3 + SC.k1) ...
                - ID.C*(ID.kd1*QDFLocal - ID.kd2*QDBLocal);
dQDBLocal = ID.C * (ID.kd1*QDFLocal - ID.kd2*QDBLocal);
dIDF = SC.k3/GC.VI(P)*QDFLocal - IDF*(ID.nDL + ID.nK) ...
           - ID.nDI*(IDF - QDF) - (SC.k1*QDF - SC.k2*QDB);
dIDB = ID.C * (ID.kd1*IDF - ID.kd2*IDB);
dQDF = -ID.nDC*QDF + ID.nDI*(IDF - QDF) - ID.C * (ID.kd1*IDF - ID.kd2*QDB);
dQDB = ID.C * (ID.kd1*QDF - ID.kd2*QDB);

% Glycaemic control model DEs.
dG = -GC.pg*(G - GFast) - SI*G*(Q+QDF)/(1 + GC.alphaG*(Q+QDF)) ... % NOTE: Removed infusion term.
          + (Ra + GC.EGP - GC.CNS)/GC.VG(P);
dI = -GC.nK(t)*I - GC.nL*I/(1 + GC.alphaI*I) - GC.nI*(I - Q) + Uen*(1 - GC.xL)/GC.VI(P); %rate of change of plasma insulin concentration
dQ = GC.nI*(I - Q) - GC.nC * Q/(1+GC.alphaG*Q);

% Pack up outgoing variables.
dY = [dqSto1;
      dqSto2;
      dqGut;
      dISC;
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