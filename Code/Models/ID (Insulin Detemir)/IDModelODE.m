function [dY] = IDModelODE(t, Y, P)
% ODE for ID model. Use with ode45.
% INPUTS:
%   t   - time at which to evaluate ODE
%   Y   - states [IDH; QDFLocal; QDBLocal; IDF; IDB; QDF; QDB] at time == t
%   P   - patient struct
% OUTPUT:
%   dY  - derivatives of states at time == t

ID = P.parameters.ID;
GC = P.parameters.GC;

%% Input
IDH       = Y(1);
QDFLocal  = Y(2);
QDBLocal  = Y(3);
IDF       = Y(4);
IDB       = Y(5);
QDF       = Y(6);
QDB       = Y(7); %NOTE: Re-removed the factor of 18...?

%% Variables
% Time dependent.
IBolus = P.data.IBolus(t); % [mU/min]

%% Computation
dIDH = -ID.ka*IDH + IBolus;
dQDFLocal = ID.ka*IDH - QDFLocal*(ID.kb + ID.kdi) ...
                - (ID.kd1*QDFLocal - ID.kd2*QDBLocal);
dQDBLocal = ID.kd1*QDFLocal - ID.kd2*QDBLocal;
dIDF = ID.kb/GC.VI*QDFLocal - IDF*(ID.nDL + ID.nK) ...
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
