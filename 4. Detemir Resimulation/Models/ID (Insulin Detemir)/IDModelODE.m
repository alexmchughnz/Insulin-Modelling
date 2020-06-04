function [dY] = IDModelODE(t, Y, P)
% ODE for ID model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [IDH; QDFLocal; QDBLocal; IDF; IDB; QDF; QDB] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

global ID

%% Input
IDH       = Y(1);
QDFLocal  = Y(2);
QDBLocal  = Y(3);
IDF       = Y(4);
IDB       = Y(5);
QDF       = Y(6);
QDB       = Y(7);

%% Variables
% Time dependent.
n = (1 + floor(t));
IBolus = P.IBolus(n);

%% Computation
dIDH = -ID.ka*IDH + IBolus;
dQDFLocal = ID.ka*IDH - QDFLocal*(ID.kb + ID.kdi) ...
                - (ID.kd1*QDFLocal - ID.kd2*QDBLocal);
dQDBLocal = ID.kd1*QDFLocal - ID.kd2*QDBLocal;
dIDF = ID.kb/GC.VI(P)*QDFLocal - IDF*(ID.nDL + ID.nK) ...
           - ID.nDI*(IDF - QDF) - (ID.kd1*QDF - ID.kd2*QDB);
dIDB = ID.kd1*IDF - ID.kd2*IDB;
dQDF = -ID.nDC*QDF + ID.nDI*(IDF - QDF) - (ID.kd1*QDF - ID.kd2*QDB);
dQDB = ID.kd1*QDF - ID.kd2*QDB;

%% Output
dY = [dIDH;
      dQDFLocal;
      dQDBLocal;
      dIDF;
      dIDB;
      dQDF;
      dQDB];

end
